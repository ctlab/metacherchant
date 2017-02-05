package tools;

import algo.MultiSequenceCalculator;
import io.RichFastaReader;
import io.graph.DeBruijnGraphUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;

import java.io.*;
import java.util.*;

public class EnvironmentFinderMultiMain extends Tool {
    public static final String NAME = "environment-finder-multi";
    public static final String DESCRIPTION = "Displays difference between multiple genomic environments";

    public static final int MAX_ENVIRONMENTS = 64;

    public final Parameter<File[]> envFiles = addParameter(new FileMVParameterBuilder("env")
            .withShortOpt("e")
            .mandatory()
            .withDescription("environment files to build difference for")
            .create());

    public final Parameter<File> seqFile = addParameter(new FileParameterBuilder("seq")
            .mandatory()
            .withDescription(".fasta file with nucleotide sequence[s]")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output")
            .mandatory()
            .withShortOpt("o")
            .withDescription("output directory to write results to")
            .create());

    private final Parameter<Integer> geneId = addParameter(new IntParameterBuilder("geneid")
            .withDefaultValue(1)
            .withShortOpt("g")
            .withDescription("gene id from .fasta file")
            .create());

    private Map<String, Integer>[] graphs;
    private String sequence;
    private String comment;
    private int k;

    private void loadInput() throws ExecutionFailedException {
        graphs = new Map[envFiles.get().length];
        for (int i = 0; i < graphs.length; i++) {
            graphs[i] = DeBruijnGraphUtils.loadGraph(envFiles.get()[i]);
        }
        if (graphs.length == 0) {
            throw new ExecutionFailedException("Zero environments given");
        }
        if (graphs.length > MAX_ENVIRONMENTS) {
            throw new ExecutionFailedException("Too many environments, expected < " + MAX_ENVIRONMENTS);
        }
        this.k = graphs[0].keySet().iterator().next().length();
        for (Map<String, Integer> graph : graphs) {
            for (String kmer : graph.keySet()) {
                if (kmer.length() != k) {
                    throw new ExecutionFailedException("K-mers of different lengths encountered: " + k + " and " + kmer.length());
                }
            }
        }

        try {
            RichFastaReader rfr = new RichFastaReader(seqFile.get());
            sequence = rfr.getDnas().get(geneId.get() - 1).toString();
            comment = rfr.getComments().get(geneId.get() - 1);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Could not load sequence file", e);
        }

    }


    protected void runImpl() throws ExecutionFailedException {
        loadInput();
        MultiSequenceCalculator calc = new MultiSequenceCalculator(sequence, k, outputDir.get().getPath(), logger, graphs);
        calc.run();
        printGene();
        logger.info("Finished processing!");
    }

    private void printGene() {
        File output = new File(outputDir.get().getPath() + "/gene.fasta");
        PrintWriter out;
        try {
            out = new PrintWriter(output);
            out.println(">" + comment);
            out.println(sequence);
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


    @Override
    protected void cleanImpl() {

    }

    public EnvironmentFinderMultiMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new EnvironmentFinderMultiMain().mainImpl(args);
    }
}
