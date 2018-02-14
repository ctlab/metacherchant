package tools;

import algo.AssemblerCalculator;
import algo.OneSequenceCalculator;
import algo.ReadsFilter;
import algo.TerminationMode;
import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import utils.FNV1AHash;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class EnvironmentAssemblerFinder extends Tool{
    public static final String NAME = "environment-assembler-finder";
    public static final String DESCRIPTION = "Finds graphic environment for many genomic sequences in given metagenomic reads in 3 stages using assembler";

    public static final int DEFAULT_MAX_THREADS = 32;

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .withDefaultValue(21)
            .create());

    public final Parameter<File[]> readsFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .withDescription("FASTQ, BINQ, FASTA reads")
            .withDefaultValue(new File[]{})
            .create());

    public final Parameter<File> seqsFile = addParameter(new FileParameterBuilder("seq")
            .mandatory()
            .withDescription("FASTA file with sequences")
            .create());


    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output")
            .mandatory()
            .withShortOpt("o")
            .withDescription("output directory")
            .create());

    public final Parameter<Integer> maxKmers = addParameter(new IntParameterBuilder("maxkmers")
//            .withDefaultValue(100000)
            .withDescription("maximum number of k-mers in created subgraph")
            .create());

    public final Parameter<Integer> maxRadius = addParameter(new IntParameterBuilder("maxradius")
//            .withDefaultValue(100000)
            .withDescription("maximum distance in k-mers from starting gene")
            .create());

    public final Parameter<Integer> minCoverage = addParameter(new IntParameterBuilder("coverage")
            .withDescription("minimum depth of k-mers to consider")
            .withDefaultValue(1)
            .create());

    public final Parameter<Boolean> bothDirections = addParameter(new BoolParameterBuilder("bothdirs")
            .withDescription("run graph search in both directions from starting sequence")
            .withDefaultValue(false)
            .create());

    public final Parameter<Integer> chunkLength = addParameter(new IntParameterBuilder("chunklength")
            .withDescription("minimum node length for BLAST search")
            .withDefaultValue(1)
            .create());

    public final Parameter<Boolean> forceHashing = addParameter(new BoolParameterBuilder("forcehash")
            .withDescription("force k-mer hashing (even for k <= 31)")
            .withDefaultValue(false)
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());

    public final Parameter<Integer> maxThreads = addParameter(new IntParameterBuilder("threads")
            .withDescription("how many java threads to use")
            .withDefaultValue(DEFAULT_MAX_THREADS)
            .create()
    );

    public final Parameter<Boolean> trimPaths = addParameter(new BoolParameterBuilder("trim")
            .withDescription("trim all not maximal paths?")
            .withDefaultValue(false)
            .create()
    );

    public final Parameter<Integer> percentFiltration = addParameter(new IntParameterBuilder("procfiltration")
            .mandatory()
            .withShortOpt("f")
            .withDescription("filtration percent // [1 .. 100]")
            .withDefaultValue(1)
            .create());

    public final Parameter<String> assembler = addParameter(new StringParameterBuilder("assembler")
            .mandatory()
            .withDescription("assembler which you want to use")
            //.withDefaultValue("None")
            .create());

    public final Parameter<String> assemblerPath = addParameter(new StringParameterBuilder("assemblerpath")
            .mandatory()
            .withDescription("path of the assembler")
            //.withDefaultValue("")
            .create());

    private BigLong2ShortHashMap reads;
    private List<DnaQ> sequences;
    private HashFunction hasher;

    public void loadInput() throws ExecutionFailedException {
        if (k.get() > 31 || forceHashing.get()) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.reads = LargeKIOUtils.loadReads(readsFiles.get(), k.get(), 0,
                    availableProcessors.get(), logger);
        } else {
            this.reads = IOUtils.loadReads(readsFiles.get(), k.get(), 0,
                    availableProcessors.get(), logger);
        }
        logger.info("Hashtable size: " + this.reads.size() + " kmers");
        try {
            this.sequences = ReadersUtils.loadDnaQs(seqsFile.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Could not load sequences from " + seqsFile.get().getPath());
        }
    }


    private HashFunction determineHashFunction() {
        if (k.get() <= 31 && !forceHashing.get()) {
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
    protected void runImpl() throws ExecutionFailedException {
        if (readsFiles.get().length > 1) {
            logger.info("EnvironmentAssemblerFinder works only with one input read!");
            return;
        }
        loadInput();

        String outputPrefix = outputDir.get().getPath() + "/";

        OneSequenceCalculator calc = new OneSequenceCalculator(sequences.get(0).toString(), k.get(),
                minCoverage.get(), outputPrefix, this.hasher, reads, logger,
                bothDirections.get(), chunkLength.get(), getTerminationMode(), trimPaths.get());
        calc.run();
        new ReadsFilter(readsFiles.get()[0], calc, outputPrefix, k.get(), percentFiltration.get(), logger).run();
        logger.info("Filtration done");
        logger.info("Finished processing all sequences!");

        new AssemblerCalculator(assembler.get(), assemblerPath.get(), outputPrefix, logger).run();
        logger.info("Finished assembling all sequences!");

        outputDir.set(new File(outputDir.get() + "/result"));
        String path = outputPrefix +
                (assembler.get().equals("spades") ? "/out_spades/contigs.fasta" : "/out_megahit/final.contigs.fasta");
        File[] f = {new File(path)};
        readsFiles.set(f);
        k.set(55);
        minCoverage.set(0);

        loadInput();

        outputPrefix = outputDir.get().getPath() + "/";

        calc = new OneSequenceCalculator(sequences.get(0).toString(), k.get(),
                minCoverage.get(), outputPrefix, this.hasher, reads, logger,
                bothDirections.get(), chunkLength.get(), getTerminationMode(), trimPaths.get());
        calc.run();
        new ReadsFilter(readsFiles.get()[0], calc, outputPrefix, k.get(), percentFiltration.get(), logger).run();
        logger.info("Filtration done");
        logger.info("Finished processing all sequences!");
    }

    @Override
    protected void cleanImpl() {

    }

    public EnvironmentAssemblerFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new EnvironmentAssemblerFinder().mainImpl(args);
    }
}
//--k 3 --reads reads.fasta --seq genes.fasta -o ress --assembler megahit --assemblerpath /home/mariya/megahit/megahit --maxkmers 1000 --coverage 0 --procfiltration 100
