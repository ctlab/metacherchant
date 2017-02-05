package io.writers;

import algo.SingleNode;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;

public class TSVWriter {
    public static final String GENE_LABEL_SUFFIX = "_start";

    private final int k;
    private final String outputPrefix;
    SingleNode[] nodes;
    private final Map<String, Integer> subgraph;
    private final int size;

    public TSVWriter(int k, String outputPrefix, SingleNode[] nodes, Map<String, Integer> subgraph) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.nodes = nodes;
        this.subgraph = subgraph;
        this.size = nodes.length;
    }

    public void print() {
        File edges = new File(outputPrefix + "/tsvs/edges.tsv");
        File nodes = new File(outputPrefix + "/tsvs/nodes.tsv");
        edges.getParentFile().mkdirs();
        printEdges(edges);
        printNodes(nodes);
    }

    private void printNodes(File outputFile) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(outputFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        out.println("id\tlength\tseq");
        for (int i = 0; i < size; i++) {
            if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id) {
                out.println((i + 1) + "\t" + nodes[i].sequence.length() + "\t" + nodes[i].sequence);
            }
        }
        out.close();
    }

    private void printEdges(File outputFile) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(outputFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        out.println("source\ttarget");
        for (SingleNode i : nodes) {
            if (!i.deleted) {
                for (SingleNode j : i.neighbors) {
                    printEdge(out, i, j);
                }
            }
        }
        out.close();
    }


    private void printEdge(PrintWriter out, SingleNode first, SingleNode second) {
        out.println(getNodeId(first) + "\t" + getNodeId(second) + "\tpp");
    }

    private String getNodeId(SingleNode node) {
        return (node.id < node.rc.id ? "" : "-") + (Math.min(node.rc.id, node.id) + 1) + (node.isGeneNode ? GENE_LABEL_SUFFIX : "");
    }
}
