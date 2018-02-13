package tools;

import algo.AssemblerCalculator;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;

import java.io.File;
import java.util.HashMap;
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
            .withDescription("assembler which you want to use")
            .withDefaultValue("None")
            .create());

    public final Parameter<String> assemblerPath = addParameter(new StringParameterBuilder("assemblerpath")
            .withDescription("path of the assembler")
            .withDefaultValue("")
            .create());


    private String getOutputPrefix(int i) {
        String outputPrefix = outputDir.get().getPath() + "/";
        if (readsFiles.get().length > 1) {
            outputPrefix += (i + 1) + "/";
        }
        return outputPrefix;
    }

    @Override
    protected void runImpl() throws ExecutionFailedException {
        /*ExecutorService execService = Executors.newFixedThreadPool(maxThreads.get());

        for (int i = 0; i < readsFiles.get().length; i++) {
            execService.execute(new AssemblerCalculator(assembler.get(), assemblerPath.get(), outputPrefix, logger));
        }
        execService.shutdown();
        */
        if (readsFiles.get().length > 1) {
            logger.info("EnvironmentAssemblerFinder works only with one input read!");
            return;
        }
        String outputPrefix = getOutputPrefix(0);
        new AssemblerCalculator(assembler.get(), assemblerPath.get(), outputPrefix, logger).run();
        logger.info("Finished assembling all sequences!");
    }

    @Override
    protected void cleanImpl() {

    }

    public EnvironmentAssemblerFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new EnvironmentFinderMain().mainImpl(args);

        Map<String, String> asmArgs = new HashMap<String, String>();
        for (int i = 0; i < args.length - 1; i += 2) {
            asmArgs.put(args[i], args[i + 1]);
        }
        String outputPrefix = "";
        if (asmArgs.containsKey("-o")) {
            asmArgs.put("-o", asmArgs.get("-o") + "/asm");
            outputPrefix = asmArgs.get("-o");
        }
        if (asmArgs.containsKey("--output")) {
            asmArgs.put("--output", asmArgs.get("--output") + "/asm");
            outputPrefix = asmArgs.get("--output");
        }
        asmArgs.put("--work-dir", asmArgs.get("--work-dir") + "/asm");
        String[] argsAssembler = new String[asmArgs.keySet().size() * 2];

        for (int i = 0; i < asmArgs.keySet().size(); i++) {
            String key = asmArgs.keySet().toArray()[i].toString();
            argsAssembler[i * 2] = key;
            argsAssembler[i * 2 + 1] = asmArgs.get(key);
        }
        new EnvironmentAssemblerFinder().mainImpl(argsAssembler);

        Map<String, String> dictionaryArgs = new HashMap<String, String>();
        for (int i = 0; i < args.length - 1; i += 2) {
            dictionaryArgs.put(args[i], args[i + 1]);
        }
        if (dictionaryArgs.get("--assembler") != null && dictionaryArgs.get("--assembler").equals("spades")) {
            dictionaryArgs.put("--reads", outputPrefix + "/out_spades/" + "contigs.fasta");
        }

        if (dictionaryArgs.get("--assembler") != null && dictionaryArgs.get("--assembler").equals("megahit")) {
            dictionaryArgs.put("--reads", outputPrefix + "/out_megahit/" + "final.contigs.fasta");
        }

        dictionaryArgs.put("--coverage", "0");
        dictionaryArgs.put("--work-dir", dictionaryArgs.get("--work-dir") + "/result");
        if (dictionaryArgs.containsKey("-o")) {
            dictionaryArgs.put("-o", dictionaryArgs.get("-o") + "/result");
        }
        if (dictionaryArgs.containsKey("--output")) {
            dictionaryArgs.put("--output", dictionaryArgs.get("--output") + "/result");
        }

        if (dictionaryArgs.containsKey("-k")) {
            dictionaryArgs.put("-k", "55");
        }
        if (dictionaryArgs.containsKey("--k")) {
            dictionaryArgs.put("--k", "55");
        }

        String[] argsAfterAssembler = new String[dictionaryArgs.keySet().size() * 2];

        for (int i = 0; i < dictionaryArgs.keySet().size(); i++) {
            String key = dictionaryArgs.keySet().toArray()[i].toString();
            argsAfterAssembler[i * 2] = key;
            argsAfterAssembler[i * 2 + 1] = dictionaryArgs.get(key);
            System.out.println(argsAfterAssembler[i] + " " + argsAfterAssembler[i+1]);
        }

        new EnvironmentFinderMain().mainImpl(argsAfterAssembler);
    }
}
//--k 3 --reads reads.fasta --seq genes.fasta -o ress --assembler megahit --assemblerpath /home/mariya/megahit/megahit --maxkmers 1000 --coverage 0 --procfiltration 100
