package io.writers;

import algo.SingleNode;
import utils.StringUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;

public class GFAWriter {
    public static final String GENE_LABEL_SUFFIX = "_start";

    private final int k;
    private final String outputPrefix;
    private SingleNode[] nodes;
    private final Map<String, Integer> subgraph;
    private final int size;
    private final String fileName;


    public GFAWriter(int k, String outputPrefix, SingleNode[] nodes, Map<String, Integer> subgraph, String fileName) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.nodes = nodes;
        this.subgraph = subgraph;
        this.size = nodes.length;
        this.fileName = fileName;
    }

    public GFAWriter(int k, String outputPrefix, SingleNode[] nodes, Map<String, Integer> subgraph) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.nodes = nodes;
        this.subgraph = subgraph;
        this.size = nodes.length;
        this.fileName = "graph";
    }

    public void print() {
        File output = new File(outputPrefix + "/" + fileName + ".gfa");
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
        for (int i = 0; i < size; i++) {
            if (!nodes[i].deleted && nodes[i].sequence.compareTo(nodes[i].rc.sequence) <= 0) {
                printLabel(out, nodes[i]);
            }
        }
        for (SingleNode i : nodes) {
            if (!i.deleted) {
                for (SingleNode j : i.neighbors) {
                    if (!j.deleted) {
                        printEdge(out, i, j);
                    }
                }
            }
        }
        out.close();
    }

    private void printEdge(PrintWriter out, SingleNode first, SingleNode second) {
        out.print("L\t");
        out.print(getNodeId(first));
        out.print('\t');
        out.print((first.sequence.compareTo(first.rc.sequence) >= 0 ? "+" : "-"));
        out.print('\t');
        out.print(getNodeId(second));
        out.print('\t');
        out.print((second.sequence.compareTo(second.rc.sequence) <= 0 ? "+" : "-"));
        out.print('\t');
        out.println((k - 1) + "M");
    }

    private String getNodeId(SingleNode node) {
        return (Math.min(node.rc.id, node.id) + 1) + (node.isGeneNode ? GENE_LABEL_SUFFIX : "");
    }

    private void printLabel(PrintWriter out, SingleNode node) {
        out.print("S\t" + getNodeId(node) + "\t" + node.sequence);
        long coverage = 0;
        for (int i = 0; i + k <= node.sequence.length(); i++) {
            String kmer = node.sequence.substring(i, i + k);
            coverage += subgraph.get(StringUtils.normalizeDna(kmer));
        }
        coverage += subgraph.get(StringUtils.normalizeDna(node.sequence.substring(node.sequence.length()-k)))*(k-1);

        out.println("\tLN:i:" + (node.sequence.length()) + "\tKC:i:" + coverage +
                (node.color == null ? "" : ("\tCL:z:" + node.color)));
    }
}
