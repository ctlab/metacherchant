package tools;

import algo.*;
import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.KmerUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;
import utils.FNV1AHash;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import static algo.SingleNode.Color.*;
import static algo.SingleNode.Color.GREY;
import static utils.StringUtils.normalizeDna;

/**
 * Created by -- on 22.01.2020.
 */
public class RecipientVisualiser extends Tool {
    public static final String NAME = "recipient-visualiser";
    public static final String DESCRIPTION = "Finds graphic environment for many genomic sequences in recipient after FMT";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> afterFiles = addParameter(new FileMVParameterBuilder("after-files")
            .mandatory()
            .withShortOpt("after")
            .withDescription("file with paired post-FMT recipient metagenomic reads")
            .create());

    public final Parameter<File> seqsFile = addParameter(new FileParameterBuilder("seq")
            .mandatory()
            .withShortOpt("seq")
            .withDescription("FASTA file with sequences")
            .create());

    public final Parameter<Integer> maxKmers = addParameter(new IntParameterBuilder("maxkmers")
//            .withDefaultValue(100000)
            .withDescription("maximum number of k-mers in created subgraph")
            .create());

    public final Parameter<Integer> maxRadius = addParameter(new IntParameterBuilder("maxradius")
            .withDefaultValue(1000)
            .withDescription("maximum distance in k-mers from starting gene")
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .optional()
            .withShortOpt("o")
            .withDescription("directory to output found reads")
            .withDefaultValue(workDir.append("graph"))
            .create());

    public final Parameter<File> inputDir = addParameter(new FileParameterBuilder("input-dir")
            .mandatory()
            .withShortOpt("i")
            .withDescription("directory containing output of reads_classifier.sh FMT classification script")
            .create());

    public final Parameter<String> extension = addParameter(new StringParameterBuilder("ext")
            .mandatory()
            .withShortOpt("ext")
            .withDescription("extension of output files of reads_classifier.sh FMT classification script")
            .create());


    private BigLong2ShortHashMap graph, from_donor, from_both, from_before, itself;
    private HashFunction hasher;
    private List<DnaQ> sequences;

    private long getKmerKey(String s) {
        if (hasher != null) {
            s = normalizeDna(s);
            return hasher.hash(s);
        } else {
            return KmerUtils.getKmerKey(DnaTools.toLong(new Dna(s)), k.get());
        }
    }

    private void loadAfterGraphs(File[] from_donor, File[] from_before, File[] from_both, File[] itself) throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.graph = LargeKIOUtils.loadReads(afterFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.from_donor = LargeKIOUtils.loadReads(from_donor, k.get(), 0, availableProcessors.get(), logger);
            this.from_before = LargeKIOUtils.loadReads(from_before, k.get(), 0, availableProcessors.get(), logger);
            this.from_both = LargeKIOUtils.loadReads(from_both, k.get(), 0, availableProcessors.get(), logger);
            this.itself = LargeKIOUtils.loadReads(itself, k.get(), 0, availableProcessors.get(), logger);
        } else {
            this.graph = IOUtils.loadReads(afterFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.from_donor = IOUtils.loadReads(from_donor, k.get(), 0, availableProcessors.get(), logger);
            this.from_before = IOUtils.loadReads(from_before, k.get(), 0, availableProcessors.get(), logger);
            this.from_both = IOUtils.loadReads(from_both, k.get(), 0, availableProcessors.get(), logger);
            this.itself = IOUtils.loadReads(itself, k.get(), 0, availableProcessors.get(), logger);
        }
        try {
            this.sequences = ReadersUtils.loadDnaQs(seqsFile.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Could not load sequences from " + afterFiles.get()[0].getPath());
        }
        //logger.info("Hashtable size: " + this.graph.size() + " kmers");
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
        sequences = null;
        from_donor = null;
        from_both = null;
        from_before = null;
        itself = null;
    }

    private Function<String, SingleNode.Color> getAfterColorNode() {
        return (seq) -> from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? RED :
                ! from_donor.contains(getKmerKey(seq)) && from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLUE :
                        ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                && from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? GREEN :
                                ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                        && ! from_both.contains(getKmerKey(seq)) && itself.contains(getKmerKey(seq)) ? YELLOW :
                                        ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                                && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLACK :
                                                GREY;
    }

    private TerminationMode getTerminationMode() throws ExecutionFailedException {
        TerminationMode termMode = new TerminationMode();
        if (maxKmers.get() == null && maxRadius.get() == null) {
            throw new ExecutionFailedException("At least one of --maxkmers and --maxradius parameters should be set");
        }
        if (maxKmers.get() != null) {
            termMode.addRestriction(TerminationMode.TerminationModeType.MAX_KMERS, maxKmers.get());
        }
        if (maxRadius.get() != null) {
            termMode.addRestriction(TerminationMode.TerminationModeType.MAX_RADIUS, maxRadius.get());
        }
        return termMode;
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        String outputPrefix = outputDir.get().getPath();
        String inputPrefix = inputDir.get().getPath() + "/";
        ExecutorService execService = Executors.newFixedThreadPool(availableProcessors.get());

        // Drawing components of after graph
        {
            logger.info("Loading after reads ...");
            File[] from_donor_files = {new File(inputPrefix + "came_from_donor_1." + extension.get()),
                    new File(inputPrefix + "came_from_donor_2." + extension.get()),
                    new File(inputPrefix + "came_from_donor_s." + extension.get())};
            File[] from_before_files = {new File(inputPrefix + "came_from_baseline_1." + extension.get()),
                    new File(inputPrefix + "came_from_baseline_2." + extension.get()),
                    new File(inputPrefix + "came_from_baseline_s." + extension.get())};
            File[] from_both_files = {new File(inputPrefix + "came_from_both_1." + extension.get()),
                    new File(inputPrefix + "came_from_both_2." + extension.get()),
                    new File(inputPrefix + "came_from_both_s." + extension.get())};
            File[] itself_files = {new File(inputPrefix + "came_itself_1." + extension.get()),
                    new File(inputPrefix + "came_itself_2." + extension.get()),
                    new File(inputPrefix + "came_itself_s." + extension.get())};
            loadAfterGraphs(from_donor_files, from_before_files, from_both_files, itself_files);

            logger.info("Creating after images ...");

            for (int i = 0; i < sequences.size(); i++) {
                execService.execute(new SeqEnvCalculator(sequences.get(i).toString(), k.get(),
                        outputPrefix + "/after", this.hasher, graph, logger,
                        getAfterColorNode(), "comp_" + i, getTerminationMode()));
            }
            execService.shutdown();
            try {
                while (!execService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS)){}
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            logger.info("Finished processing all sequences!");
        }
    }

    public RecipientVisualiser() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new RecipientVisualiser().mainImpl(args);
    }
}