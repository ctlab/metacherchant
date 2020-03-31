package tools;

import algo.TripleFinder;
import algo.TripleFinder2;
import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.io.sources.PairSource;
import ru.ifmo.genetics.io.writers.WritersUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.tools.io.LazyDnaQReaderTool;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.pairs.UniPair;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import utils.FNV1AHash;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Created by -- on 18.03.2019.
 */
public class TripleReadsClassifier extends Tool {
    public static final String NAME = "triple-reads-classifier";
    public static final String DESCRIPTION = "classifies reads based on weighted De Bruijn graph with two values of k-mers" +
            " and splits them into three categories";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<Integer> k2 = addParameter(new IntParameterBuilder("k2")
            .mandatory()
            .withShortOpt("k2")
            .withDescription("second k-mer size. k2 > k")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("input-files")
            .withShortOpt("i")
            .withDescription("file with paired input reads for De Bruijn graph")
            .create());

    public final Parameter<File[]> inputKmers1 = addParameter(new FileMVParameterBuilder("input-kmers-1")
            .withShortOpt("ik1")
            .withDescription("file with k-mers in binary format for De Bruijn graph")
            .create());

    public final Parameter<File[]> inputKmers2 = addParameter(new FileMVParameterBuilder("input-kmers-2")
            .withShortOpt("ik2")
            .withDescription("file with k-mers in binary format for De Bruijn graph")
            .create());

    public final Parameter<File[]> readsFiles = addParameter(new FileMVParameterBuilder("read-files")
            .mandatory()
            .withShortOpt("r")
            .withDescription("files with paired reads to classify")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .optional()
            .withShortOpt("o")
            .withDescription("directory to output found reads")
            .withDefaultValue(workDir.append("reads_classifier"))
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());

    public final Parameter<Boolean> doCorrection = addParameter(new BoolParameterBuilder("correction")
            .optional()
            .withShortOpt("corr")
            .withDescription("Do replacement of nucleotide in read with one low quality position")
            .withDefaultValue(false)
            .create());

    public final Parameter<Boolean> interval95 = addParameter(new BoolParameterBuilder("interval95")
            .optional()
            .withDescription("Set the interval width to probability 0.95")
            .withDefaultValue(false)
            .create());

    public final Parameter<Integer> found_threshold = addParameter(new IntParameterBuilder("found-threshold")
            .optional()
            .withShortOpt("found")
            .withDescription("Minimum coverage breadth for class `found` [0 - 100 %]")
            .withDefaultValue(90)
            .create());

    public final Parameter<Integer> half_threshold = addParameter(new IntParameterBuilder("half-threshold")
            .optional()
            .withShortOpt("half")
            .withDescription("Minimum coverage breadth for class `half-found` [0 - 100 %]")
            .withDefaultValue(40)
            .create());


    private BigLong2ShortHashMap graph;
    private HashFunction hasher;

    public TripleReadsClassifier() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new TripleReadsClassifier().mainImpl(args);
    }

    @Override
    protected void cleanImpl() {
        graph = null;
        hasher = null;
    }

    public void loadGraph(int k, Parameter<File[]> inputKmers) throws ExecutionFailedException {
        if (k > 31) {
            this.hasher = LargeKIOUtils.hash = determineHashFunction(k);
        }
        if (inputKmers.get() != null && inputKmers.get()[0].getName().toLowerCase().endsWith("kmers.bin")) {
            this.graph = IOUtils.loadKmers(inputKmers.get(), 0, availableProcessors.get(), logger);
        }
        else {
            if (k > 31) {
                logger.info("Reading hashes of k-mers instead");
                this.graph = LargeKIOUtils.loadReads(inputFiles.get(), k, 0,
                        availableProcessors.get(), logger);
            } else {
                this.graph = IOUtils.loadReads(inputFiles.get(), k, 0,
                        availableProcessors.get(), logger);
            }
        }
        logger.info("Hashtable size: " + this.graph.size() + " kmers");
    }

    private HashFunction determineHashFunction(int k) {
        if (k <= 31) {
            return null;
        }
        String name = hashFunction.get().toLowerCase();
        if (name.equals("fnv1a")) {
            logger.info("Using FNV1a hash function");
            return new FNV1AHash();
        } else {
            logger.info("Using default polynomial hash function");
            return new PolynomialHash();
        }
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        if (k.get() >= k2.get()) {
            throw new ExecutionFailedException("k2 should be greater than k, given: " + k.get() + " " + k2.get());
        }
        outputDir.get().mkdirs();

        info("Loading reads...");
        LazyDnaQReaderTool dnaQReader = new LazyDnaQReaderTool();
        ArrayList<NamedSource<? extends LightDnaQ>> sources = new ArrayList<>(readsFiles.get().length);
        for (File file : readsFiles.get()) {
            dnaQReader.fileIn.set(file);
            dnaQReader.simpleRun();
            NamedSource<? extends LightDnaQ> source = dnaQReader.dnaQsSourceOut.get();
            sources.add(source);
        }
        NamedSource<? extends LightDnaQ> source1 = sources.get(0);
        NamedSource<? extends LightDnaQ> source2 = sources.get(1);
        PairSource<LightDnaQ> pairedSource = PairSource.create(source1, source2);

        info("Building graph with k = " + k.get() + " ...");
        loadGraph(k.get(), inputKmers1);
        info(doCorrection.get() ? "Searching for corrected reads in graph..." : "Searching for reads in graph...");

        Map<String, FindResult> isFoundInGraphOne_1 = new ConcurrentHashMap<>();
        Map<String, FindResult> isFoundInGraphOne_2 = new ConcurrentHashMap<>();
        ExecutorService executorService = Executors.newFixedThreadPool(availableProcessors.get());
        for (UniPair<LightDnaQ> pair : pairedSource) {
            executorService.execute(new TripleFinder(pair, k.get(), graph, hasher, isFoundInGraphOne_1, isFoundInGraphOne_2,
                    doCorrection.get(), interval95.get() ? 1.96 : 1,(double)found_threshold.get() / 100,
                    (double)half_threshold.get() / 100));
        }
        executorService.shutdown();
        try {
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            throw new ExecutionFailedException("Error while searching reads in graph: " + e.toString());
        }

        cleanImpl();
        info("Building graph with k = " + k2.get() + " ...");
        loadGraph(k2.get(), inputKmers2);
        info(doCorrection.get() ? "Searching for corrected reads in graph..." : "Searching for reads in graph...");

        Queue<UniPair<LightDnaQ>> both_found = new ConcurrentLinkedQueue<>();
        Queue<UniPair<LightDnaQ>> both_half_found = new ConcurrentLinkedQueue<>();
        Queue<UniPair<LightDnaQ>> both_not_found = new ConcurrentLinkedQueue<>();
        Queue<LightDnaQ> s_found = new ConcurrentLinkedQueue<>();
        Queue<LightDnaQ> s_half_found = new ConcurrentLinkedQueue<>();
        Queue<LightDnaQ> s_not_found = new ConcurrentLinkedQueue<>();

        executorService = Executors.newFixedThreadPool(availableProcessors.get());
        for (UniPair<LightDnaQ> pair : pairedSource) {
            executorService.execute(new TripleFinder2(pair, k2.get(), graph, hasher,
                    isFoundInGraphOne_1, isFoundInGraphOne_2, doCorrection.get(),
                    both_found, both_half_found, both_not_found, s_found, s_half_found, s_not_found,
                    interval95.get() ? 1.96 : 1, (double)found_threshold.get() / 100,
                    (double)half_threshold.get() / 100));
        }
        executorService.shutdown();
        try {
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            throw new ExecutionFailedException("Error while searching reads in graph: " + e.toString());
        }

        FoundStats stats = new FoundStats(both_found.size(), both_half_found.size(), both_not_found.size(),
                s_found.size(), s_half_found.size(), s_not_found.size());

        debug("Free memory left: " + Misc.availableMemoryWithoutRunningGCAsString());
        info("|\t" + "Total: " + stats.getTotal() + " reads");
        info("|\t" + "Paired: " + stats.getPaired() + " reads");
        info("|\t" + "Total quality: " + String.format("%.2f", 100 * (double) stats.getPaired() / stats.getTotal()) + " %");
        info("|\t" + "Found: " + stats.getFound() + " reads");
        info("|\t" + "Percent of found reads: " + String.format("%.2f", 100 * (double) stats.getFound() / stats.getTotal()) + " %");
        info("|\t" + "Quality of found bin: " + String.format("%.2f", stats.getQualityFound()) + " %");
        info("|\t" + "Not found: " + stats.getNotFound() + " reads");
        info("|\t" + "Percent of not found reads: " + String.format("%.2f", 100 * (double) stats.getNotFound() / stats.getTotal()) + " %");
        info("|\t" + "Quality of not found bin: " + String.format("%.2f", stats.getQualityNotFound()) + " %");

        info("|\t" + "Half found: " + stats.getHalfFound() + " reads");
        info("|\t" + "Percent of half found reads: " + String.format("%.2f", 100 * (double) stats.getHalfFound() / stats.getTotal()) + " %");
        info("|\t" + "Quality of half found bin: " + String.format("%.2f", stats.getQualityHalfFound()) + " %");


        info("Writing classified reads...");
        WritersUtils.writeDnaQsToFastqFile(both_found.stream().map(UniPair::first)
                .collect(Collectors.toList()), new File(outputDir.get() + "/found_1.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_found.stream().map(UniPair::second)
                .collect(Collectors.toList()), new File(outputDir.get() + "/found_2.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_half_found.stream().map(UniPair::first).filter(a -> a.length() > 0)
                .collect(Collectors.toList()), new File(outputDir.get() + "/half_found_1.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_half_found.stream().map(UniPair::second).filter(a -> a.length() > 0)
                .collect(Collectors.toList()), new File(outputDir.get() + "/half_found_2.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_not_found.stream().map(UniPair::first).filter(a -> a.length() > 0)
                .collect(Collectors.toList()), new File(outputDir.get() + "/not_found_1.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_not_found.stream().map(UniPair::second).filter(a -> a.length() > 0)
                .collect(Collectors.toList()), new File(outputDir.get() + "/not_found_2.fastq"));

        WritersUtils.writeDnaQsToFastqFile(s_found.stream().filter(a -> a.length() > 0).collect(Collectors.toList()),
                new File(outputDir.get() + "/found_s.fastq"));
        WritersUtils.writeDnaQsToFastqFile(s_half_found.stream().filter(a -> a.length() > 0).collect(Collectors.toList()),
                new File(outputDir.get() + "/half_found_s.fastq"));
        WritersUtils.writeDnaQsToFastqFile(s_not_found.stream().filter(a -> a.length() > 0).collect(Collectors.toList()),
                new File(outputDir.get() + "/not_found_s.fastq"));

        info("Reads have been written. Finishing...");
    }

    public enum FindResult {
        FOUND, HALF_FOUND, NOT_FOUND
    }

    class FoundStats {
        private int both_found;
        private int both_notFound;
        private int both_halfFound;
        private int found;
        private int notFound;
        private int halfFound;

        FoundStats() {
            both_found = 0;
            both_notFound = 0;
            both_halfFound = 0;
            found = 0;
            notFound = 0;
            halfFound = 0;
        }

        FoundStats(int both_found, int both_halfFound, int both_notFound, int found, int halfFound, int notFound) {
            this.both_found = both_found;
            this.both_halfFound = both_halfFound;
            this.both_notFound = both_notFound;
            this.found = found;
            this.halfFound = halfFound;
            this.notFound = notFound;
        }

        int getTotal() {
            return 2 * (both_notFound + both_found + both_halfFound) + found + notFound + halfFound;
        }

        int getHalfFound() {
            return 2 * both_halfFound + halfFound;
        }

        int getFound() {
            return 2 * both_found + found;
        }

        int getNotFound() {
            return 2 * both_notFound + notFound;
        }

        int getPaired() {
            return 2 * (both_found + both_notFound + both_halfFound);
        }

        double getQualityFound() {
            return (double) both_found * 2 / getFound() * 100;
        }

        double getQualityNotFound() {
            return (double) both_notFound * 2 / getNotFound() * 100;
        }

        double getQualityHalfFound() {
            return (double) both_halfFound * 2 / getHalfFound() * 100;
        }
    }

}
