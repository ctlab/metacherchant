package io.writers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Set;

public class GraphvizWriter {
    public static final String GENE_LABEL_SUFFIX = "_start";

    private final int k;
    private final String outputPrefix;
    private final boolean[] deleted;
    private final String[] label;
    private final boolean[] isGeneNode;
    private final Map<String, Integer> subgraph;
    private final Set<Integer>[] inEdges;
    private final Set<Integer>[] outEdges;
    private final int size;


    public GraphvizWriter(int k, String outputPrefix, boolean[] deleted, String[] label, boolean[] isGeneNode, Map<String, Integer> subgraph,
                          Set<Integer>[] inEdges, Set<Integer>[] outEdges) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.deleted = deleted;
        this.label = label;
        this.isGeneNode = isGeneNode;
        this.subgraph = subgraph;
        this.inEdges = inEdges;
        this.outEdges = outEdges;
        this.size = inEdges.length;
    }

    public void print() {
        File outputDot = new File(outputPrefix + "/graph.dot");
        File outputPng = new File(outputPrefix + "/graph.png");
        outputDot.getParentFile().mkdirs();
        printGraph(outputDot);
        makePng(outputDot, outputPng);
    }


    private void printGraph(File outputFile) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(outputFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        out.println("strict digraph G {");
        for (int i = 0; i < size; i++) {
            if (!deleted[i]) {
                printLabel(out, label[i], i);
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j : outEdges[i]) {
                if (!deleted[i] && !deleted[j]) {
                    printEdge(out, i + 1, j + 1, label[j].charAt(k - 1));
                }
            }
        }
        out.println("}");
        out.close();
    }

    private void printEdge(PrintWriter out, int id1, int id2, char label) {
        out.println(id1 + "->" + id2 + " [label=\" " + label + "\"]");
    }

    private void printLabel(PrintWriter out, String label, int id) {
        out.print((id + 1) + " [");
        if (isGeneNode[id]) {
            out.print("fontcolor=red style=\"bold\" label=<<B>");
            out.print(label.length());
            printFreqs(out, label);
            out.print("</B>>");
            out.println("]");
            return;
        }
        out.print("label=<");
        out.print(label.length());
        printFreqs(out, label);
        out.print(">");
        out.println("]");
    }

    private void printFreqs(PrintWriter out, String label) {
        int min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
        for (int i = 0; i + k <= label.length(); i++) {
            String kmer = label.substring(i, i + k);
            int freq = subgraph.get(kmer);
            min = Math.min(min, freq);
            max = Math.max(max, freq);
        }
        out.printf(":%d:%d", min, max);
    }

    private void makePng(File graph, File png) {
        try {
            Process p = Runtime.getRuntime().exec(new String[]{"dot", "-Tpng", graph.getPath(), "-o", png.getPath()});
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
