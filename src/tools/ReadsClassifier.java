package tools;

import algo.PairFinder;
import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.io.sources.PairSource;
import ru.ifmo.genetics.io.writers.WritersUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.tools.io.LazyDnaQReaderTool;
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
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * Created by -- on 11.03.2019.
 */
public class ReadsClassifier extends Tool {
    public static final String NAME = "reads-classifier";
    public static final String DESCRIPTION = "classifies reads based on weighted De Bruijn graph";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            //.withDefaultValue(21)
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("input-files")
            .mandatory()
            .withShortOpt("i")
            .withDescription("file with paired input reads for De Bruijn graph")
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

    public final Parameter<Boolean> doCorrection = addParameter(new BoolParameterBuilder("correction")
            .optional()
            .withShortOpt("corr")
            .withDescription("Do replacement of nucleotide in read with one low quality position")
            .withDefaultValue(false)
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
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


    private BigLong2ShortHashMap graph;
    private HashFunction hasher;

    public void loadGraph() throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.graph = LargeKIOUtils.loadReads(inputFiles.get(), k.get(), 0,
                    availableProcessors.get(), logger);
        } else {
            this.graph = IOUtils.loadReads(inputFiles.get(), k.get(), 0,
                    availableProcessors.get(), logger);
        }
        logger.info("Hashtable size: " + this.graph.size() + " kmers");
    }


    private HashFunction determineHashFunction() {
        if (k.get() <= 31) {
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
    protected void cleanImpl() {
        graph = null;
        hasher = null;
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        outputDir.get().mkdirs();
        loadGraph();

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

        Queue<UniPair<LightDnaQ>> both_found = new ConcurrentLinkedQueue<>();
        Queue<UniPair<LightDnaQ>> first_found = new ConcurrentLinkedQueue<>();
        Queue<UniPair<LightDnaQ>> second_found = new ConcurrentLinkedQueue<>();
        Queue<UniPair<LightDnaQ>> both_not_found = new ConcurrentLinkedQueue<>();

        info(doCorrection.get() ? "Searching for corrected reads in graph..." : "Searching for reads in graph...");

        ExecutorService executorService = Executors.newFixedThreadPool(availableProcessors.get());
        for (UniPair<LightDnaQ> pair : pairedSource) {
           executorService.execute(new PairFinder(pair, k.get(), graph, hasher, both_found, first_found,
                   second_found, both_not_found, doCorrection.get(), interval95.get() ? 1.96 : 1,
                   (double)found_threshold.get() / 100));
        }
        executorService.shutdown();
        try {
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            throw new ExecutionFailedException("Error while searching reads in graph: " + e.toString());
        }

        FoundStats stats = new FoundStats(both_found.size(), first_found.size(), second_found.size(), both_not_found.size());

        info("|\t" + "Total: " + stats.getTotal() + " reads");
        info("|\t" + "Paired: " + stats.getPaired() + " reads");
        info("|\t" + "Total quality: " + String.format("%.2f", 100 * (double)stats.getPaired() / stats.getTotal()) + " %");
        info("|\t" + "Found: " + stats.getFound() + " reads");
        info("|\t" + "Percent of found reads: " + String.format("%.2f", 100 * (double)stats.getFound() / stats.getTotal()) + " %");
        info("|\t" + "Quality of found bin: " + String.format("%.2f", stats.getQualityFound()) + " %");
        info("|\t" + "Not found: " + stats.getNotFound() + " reads");
        info("|\t" + "Percent of not found reads: " + String.format("%.2f", 100 * (double)stats.getNotFound() / stats.getTotal()) + " %");
        info("|\t" + "Quality of not found bin: " + String.format("%.2f", stats.getQualityNotFound()) + " %");


        /*if (stats.getTotal() != found_reads_1.size() + found_reads_2.size() + not_found_reads_1.size() +
                not_found_reads_2.size() + found_reads_s.size() + not_found_reads_s.size()) {
            throw new ExecutionFailedException("Total counted badly!");
        }*/
        info("Writing classified reads...");
        WritersUtils.writeDnaQsToFastqFile(both_found.stream().map(UniPair::first).collect(Collectors.toList()),
                                            new File(outputDir.get() + "/found_1.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_found.stream().map(UniPair::second).collect(Collectors.toList()),
                                            new File(outputDir.get() + "/found_2.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_not_found.stream().map(UniPair::first).collect(Collectors.toList()),
                                            new File(outputDir.get() + "/not_found_1.fastq"));
        WritersUtils.writeDnaQsToFastqFile(both_not_found.stream().map(UniPair::second).collect(Collectors.toList()),
                                            new File(outputDir.get() + "/not_found_2.fastq"));

        List<LightDnaQ> found_s = first_found.stream().map(UniPair::first).filter(a -> a.length() > 0).collect(Collectors.toList());
        found_s.addAll(second_found.stream().map(UniPair::second).filter(a -> a.length() > 0).collect(Collectors.toList()));
        WritersUtils.writeDnaQsToFastqFile(found_s, new File(outputDir.get() + "/found_s.fastq"));
        List<LightDnaQ> not_found_s = first_found.stream().map(UniPair::second).filter(a -> a.length() > 0).collect(Collectors.toList());
        not_found_s.addAll(second_found.stream().map(UniPair::first).filter(a -> a.length() > 0).collect(Collectors.toList()));
        WritersUtils.writeDnaQsToFastqFile(not_found_s, new File(outputDir.get() + "/not_found_s.fastq"));
        info("Reads have been written. Finishing...");
    }

    class FoundStats {
        private int both_found;
        private int first_found_second_notFound;
        private int first_notFound_second_found;
        private int both_notFound;

        FoundStats() {
            both_found = 0;
            first_found_second_notFound = 0;
            first_notFound_second_found = 0;
            both_notFound = 0;
        }

        FoundStats(int both_found, int first_found_second_notFound, int first_notFound_second_found, int both_notFound) {
            this.both_found = both_found;
            this.first_found_second_notFound = first_found_second_notFound;
            this.first_notFound_second_found = first_notFound_second_found;
            this.both_notFound = both_notFound;
        }

        int getTotal() {
            return 2 * (both_notFound + both_found + first_found_second_notFound + first_notFound_second_found);
        }

        int getFound() {
            return 2 * both_found + first_found_second_notFound + first_notFound_second_found;
        }

        int getNotFound() {
            return 2 * both_notFound + first_found_second_notFound + first_notFound_second_found;
        }

        int getPaired() {
            return 2 * (both_found + both_notFound);
        }

        double getQualityFound() {
            return (double)both_found * 2 / (both_found * 2 + first_found_second_notFound + first_notFound_second_found) * 100;
        }

        double getQualityNotFound() {
            return (double)both_notFound * 2 / (both_notFound * 2 + first_found_second_notFound + first_notFound_second_found) * 100;
        }
    }

    public ReadsClassifier() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new ReadsClassifier().mainImpl(args);
    }
}
