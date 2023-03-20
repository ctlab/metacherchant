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

    public static final int MAX_ENVIRONMENTS = 256;

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
            logger.warn("Found more than 256 environments. Grayscale graph may be not accurate.");
            //throw new ExecutionFailedException("Too many environments, expected < " + MAX_ENVIRONMENTS);
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
        printProbability();
        
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
    
    private void printProbability(){
    	 File outputSym = new File(outputDir.get().getPath() + "/Jacard_sym.txt");
    	 File outputAlt = new File(outputDir.get().getPath() + "/Jacard_alt.txt");
         PrintWriter outSym;
         PrintWriter outAlt;
         try {
             outSym = new PrintWriter(outputSym);
             outAlt = new PrintWriter(outputAlt);
             
             outSym.println("The" + "[31mWarning! " + "symmetric <<Jaccard distance>> (1 - AB/AUB):");
             outSym.println();
             outAlt.println("The" + "[31mWarning! " + "alternative <<Jaccard distance>> (1 - AB/A):");
             outAlt.println();
        	 int[][] difference = new int[graphs.length][graphs.length];
        	 int[][] differenceAlt = new int[graphs.length][graphs.length];
             int[][] union = new int[graphs.length][graphs.length]; 
             int i = 0;
             for (Map<String, Integer> graphF: graphs){
             	int j = 0;
             	for (Map<String, Integer> graphS: graphs){
             		for (String kmer : graphF.keySet()) { 
             			if (graphS.containsKey(kmer) == false){
             				difference[i][j] += graphF.get(kmer);
             				differenceAlt[i][j] += graphF.get(kmer);
             				union[i][j] += graphF.get(kmer);  
             			}
             			else{
             				difference[i][j] += Math.abs(graphF.get(kmer) - graphS.get(kmer));
             				differenceAlt[i][j] += Math.abs(graphF.get(kmer) - graphS.get(kmer));
             				union[i][j] += Math.max(graphF.get(kmer), graphS.get(kmer));
             			}
             				
             		}
             		for (String kmer : graphS.keySet()) { 
             			if (graphF.containsKey(kmer) == false){
             				difference[i][j] += graphS.get(kmer);
             				union[i][j] += graphS.get(kmer);  
             			}
             		}
             		j++;
             	}
             	i++;
             }
             
             for (i = 0; i < graphs.length; i++){
             	outSym.print(envFiles.get()[i]);
             	outAlt.print(envFiles.get()[i]);
             	for (int j = 0; j < graphs.length; j++){
             		int intersection = union[i][j] - difference[i][j];
             		outSym.printf("%6.2f ", 1 - intersection/(float)union[i][j]);
             		outAlt.printf("%6.2f ", 1 - intersection/(float)(union[i][j] - differenceAlt[i][j]));
             	}
             	outSym.println();
             	outAlt.println();
             }
             
             outSym.close();
             outAlt.close();
        	 
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
