package io.writers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Set;

public class GMLWriter {
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


    public GMLWriter(int k, String outputPrefix, boolean[] deleted, String[] label, boolean[] isGeneNode, Map<String, Integer> subgraph,
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
        File output = new File(outputPrefix + "/graph.gml");
        output.getParentFile().mkdirs();
        printGraph(output);
    }


    private void printGraph(File outputFile) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(outputFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        out.println("writers [");
        for (int i = 0; i < size; i++) {
            if (!deleted[i]) {
                printNode(out, label[i], i);
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j : outEdges[i]) {
                if (!deleted[i] && !deleted[j]) {
                    printEdge(out, i, j, label[j].charAt(k - 1));
                }
            }
        }
        out.println("]");
        out.close();
    }

    private void printEdge(PrintWriter out, int src, int dest, char label) {
        out.println("edge [");
        out.println("source \"" + getNodeId(src) + '"');
        out.println("target \"" + getNodeId(dest) + '"');
        out.println("]");
    }

    private String getNodeId(int id) {
        return (id + 1) + (isGeneNode[id] ? GENE_LABEL_SUFFIX : "");
    }

    private void printNode(PrintWriter out, String label, int id) {
        out.println("node [");

        out.println("id \"" + getNodeId(id) + '"');
        out.println("length " + label.length());
        out.println("seq \"" + label + '"');

        out.println("]");
    }
}
