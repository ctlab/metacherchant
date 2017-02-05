package io.writers;

import algo.SingleNode;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;

public class LastGraphWriter {
    private final int k;
    private final String outputPrefix;
    private final SingleNode[] nodes;
    private final Map<String, Integer> subgraph;
    private final int size;


    public LastGraphWriter(int k, String outputPrefix, SingleNode[] nodes, Map<String, Integer> subgraph) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.nodes = nodes;
        this.subgraph = subgraph;
        this.size = nodes.length;
    }

    public void print() {
        File output = new File(outputPrefix + "graph_LastGraph");
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
        int realNodes = 0;
        for (int i = 0; i < size; i++) {
            if (!nodes[i].deleted) {
                realNodes++;
            }
        }
        out.println(realNodes + " 0 " + k + " " + 1);
        for (int i = 0; i < size; i++) {
            if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id) {
                printLabel(out, nodes[i]);
            }
        }
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
        out.print("ARC ");
        out.print((first.id < first.rc.id ? "" : "-") + (Math.min(first.rc.id, first.id) + 1) + (first.isGeneNode ? "_oxa347" : ""));
        out.print(' ');
        out.print((second.id > second.rc.id ? "" : "-") + (Math.min(second.rc.id, second.id) + 1) + (second.isGeneNode ? "_oxa347" : ""));
        out.println();
    }

    private String getNodeId(SingleNode node) {
        return (node.id < node.rc.id ? "" : "-") + (Math.min(node.rc.id, node.id) + 1) + (node.isGeneNode ? "_oxa347" : "");
    }

    private void printLabel(PrintWriter out, SingleNode node) {
        out.print("NODE " + getNodeId(node));
        long coverage = 0;
        for (int i = 0; i + k <= node.sequence.length(); i++) {
            String kmer = node.sequence.substring(i, i + k);
            coverage += subgraph.get(kmer);
        }
        out.println(" " + (node.sequence.length()) + " " + coverage + " " + coverage + " 0 0");

        out.println(node.sequence);
        out.println(node.rc.sequence);
    }
}
