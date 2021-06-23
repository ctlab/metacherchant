package tools;

import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
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

/**
 * Created by -- on 31.03.2020.
 */
public class KmersCounter extends Tool {
    public static final String NAME = "kmer-counter";
    public static final String DESCRIPTION = "Count k-mers in given reads with BigLong2ShortHashMap";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of reads files from single environment. FASTQ, BINQ, FASTA (ignored reads with 'N')")
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());



    private BigLong2ShortHashMap graph;

    public void loadGraph() throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            LargeKIOUtils.hash = determineHashFunction();
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
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        outputDir.get().mkdirs();
        loadGraph();
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        String name = ReadersUtils.readDnaLazy(inputFiles.get()[0]).name();
        File outFile = new File(outputDir.get(), name + ".kmers.bin");
        File stFile = new File(outputDir.get(), name + ".stat.txt");


        debug("Starting to print k-mers to " + outFile.getPath());
        long c = 0;
        try {
            c = IOUtils.printKmers(graph, 0, outFile, stFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        info(NumUtils.groupDigits(graph.size()) + " k-mers found, "
                + NumUtils.groupDigits(c) + " (" + String.format("%.1f", c * 100.0 / graph.size()) + "%) of them is good (not erroneous)");

        if (graph.size() == 0) {
            warn("No k-mers found in reads! Perhaps you reads file is empty or k-mer size is too big");
        } else if (c == 0 || c < (long) (graph.size() * 0.03)) {
            warn("Too few good k-mers were found! Perhaps you should decrease k-mer size or --maximal-bad-frequency value");
        }
        long allKmersNumber = (1L << (2*k.get())) / 2;  // (4^k)/2
        if (graph.size() == allKmersNumber) {
            warn("All possible k-mers were found in reads! Perhaps you should increase k-mer size");
        } else if (graph.size() >= (long) (allKmersNumber * 0.99)) {
            warn("Almost all possible k-mers were found in reads! Perhaps you should increase k-mer size");
        }

        info("k-mers printed to " + outFile.getPath());
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new KmersCounter().mainImpl(args);
    }

    public KmersCounter() {
        super(NAME, DESCRIPTION);
    }

}
