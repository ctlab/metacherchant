package tools;

import io.IOUtils;
import io.LargeKIOUtils;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.tools.io.LazyDnaQReaderTool;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import utils.FNV1AHash;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created by -- on 06.04.2020.
 */
public class SequenceCoverage extends Tool {
    public static final String NAME = "seq-cov";
    public static final String DESCRIPTION = "Calculates coverage of sequences by k-mers from metagenomic bins";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> beforeFiles = addParameter(new FileMVParameterBuilder("from-before")
            .mandatory()
            .withDescription("file with paired input reads for came_from_before bin")
            .create());

    public final Parameter<File[]> donorFiles = addParameter(new FileMVParameterBuilder("from-donor")
            .mandatory()
            .withDescription("file with paired input reads for came_from_donor bin")
            .create());

    public final Parameter<File[]> bothFiles = addParameter(new FileMVParameterBuilder("from-both")
            .mandatory()
            .withDescription("file with paired input reads for came_from_both bin")
            .create());

    public final Parameter<File[]> itselfFiles = addParameter(new FileMVParameterBuilder("itself")
            .mandatory()
            .withDescription("file with paired input reads for came_itself bin")
            .create());

    public final Parameter<File> seqFile = addParameter(new FileParameterBuilder("read-file")
            .mandatory()
            .withShortOpt("r")
            .withDescription("file with sequences to classify")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .optional()
            .withShortOpt("o")
            .withDescription("directory to output found reads")
            .withDefaultValue(workDir.append("sequence_coverage"))
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());


    private BigLong2ShortHashMap donor, before, both, itself;
    private HashFunction hasher;

    public SequenceCoverage() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new SequenceCoverage().mainImpl(args);
    }

    public BigLong2ShortHashMap loadGraph(File[] inputFiles) throws ExecutionFailedException {
        BigLong2ShortHashMap graph;
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            graph = LargeKIOUtils.loadReads(inputFiles, k.get(), 0,
                    availableProcessors.get(), logger);
        } else {
            graph = IOUtils.loadReads(inputFiles, k.get(), 0,
                    availableProcessors.get(), logger);
        }
        logger.info("Hashtable size: " + graph.size() + " kmers");
        return graph;
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
        donor = null;
        before = null;
        both = null;
        itself = null;
        hasher = null;
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        outputDir.get().mkdirs();

        info("Loading bins ...");
        donor = loadGraph(donorFiles.get());
        before = loadGraph(beforeFiles.get());
        both = loadGraph(bothFiles.get());
        itself = loadGraph(itselfFiles.get());

        info("Calculating sequence coverage...");
        LazyDnaQReaderTool dnaQReader = new LazyDnaQReaderTool();
        dnaQReader.fileIn.set(seqFile.get());
        dnaQReader.simpleRun();
        NamedSource<? extends LightDnaQ> source = dnaQReader.dnaQsSourceOut.get();

        int k = this.k.get();
        FileWriter out = new FileWriter(outputDir.get() + "/seq_cov.csv");
        PrintWriter printer = new PrintWriter(out);
        printer.println("name, from_donor_depth, from_donor_breadth, from_before_depth, from_before_breadth" +
                ", from_both_depth, from_both_breadth, itself_depth, itself_breadth");

        for (LightDnaQ seq : source) {
            Dna dna = new Dna(seq);
            printer.print(dna.toString());
            printSeqBin(donor, dna, k, printer);
            printSeqBin(before, dna, k, printer);
            printSeqBin(both, dna, k, printer);
            printSeqBin(itself, dna, k, printer);
            printer.println();
        }

        printer.close();
        info("Processed all sequences...");

    }

    private void printSeqBin(BigLong2ShortHashMap graph, Dna seq, int k, PrintWriter printer) {
        long depth = 0, breadth = 0;

        if (k > 31) {
            for (int i = 0; i + k <= seq.length(); i++) {
                long hash = hasher.hash(seq, i, i + k);
                short value = graph.getWithZero(hash);
                if (value > 0) {
                    breadth++;
                }
                depth += value;
            }
        } else {
            for (ShortKmer kmer : ShortKmer.kmersOf(seq, k)) {
                short value = graph.getWithZero(kmer.toLong());
                if (value > 0) {
                    breadth++;
                }
                depth += value;
            }
        }

        printer.print(", " + depth * 1. / (seq.length() - k + 1) + ", " + breadth * 1. / (seq.length() - k + 1));
    }
}
