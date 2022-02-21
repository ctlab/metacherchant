package io.writers;

import algo.MultiNode;
import utils.StringUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;

public class GFAWriterMulti {
    public static final String GENE_LABEL_SUFFIX = "_start";

    public static final String COLOR_BLACK = "#000000";
    public static final String COLOR_RED = "#ff0000";
    public static final String COLOR_GREEN = "#00ff00";
    public static final String COLOR_BLUE = "#0000ff";

    private final int k;
    private final String outputPrefix;
    private final MultiNode[] nodes;
    private final int size;
    private final Map<String, Integer>[] graphs;

    public GFAWriterMulti(int k, String outputPrefix, MultiNode[] nodes, Map<String, Integer>[] graphs) {
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.nodes = nodes;
        this.size = nodes.length;
        this.graphs = graphs;
    }

    public void print() {
        File output = new File(outputPrefix + "/graph.gfa");
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
            if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id) {
                printLabel(out, nodes[i]);
            }
        }
        for (MultiNode i : nodes) {
            if (!i.deleted) {
                for (MultiNode j : i.neighbors) {
                    printEdge(out, i, j);
                }
            }
        }
        out.close();
    }

    private void printEdge(PrintWriter out, MultiNode first, MultiNode second) {
        out.print("L\t");
        out.print((Math.min(first.rc.id, first.id) + 1) + (first.isGeneNode ? GENE_LABEL_SUFFIX : ""));
        out.print('\t');
        out.print((first.id < first.rc.id ? "+" : "-"));
        out.print('\t');
        out.print((Math.min(second.rc.id, second.id) + 1) + (second.isGeneNode ? GENE_LABEL_SUFFIX : ""));
        out.print('\t');
        out.print((second.id > second.rc.id ? "+" : "-"));
        out.print('\t');
        out.println((k - 1) + "M");
    }

    private String getNodeId(MultiNode node) {
        return (node.id < node.rc.id ? "" : "-") + (Math.min(node.rc.id, node.id) + 1) + (node.isGeneNode ? GENE_LABEL_SUFFIX : "");
    }

    private void printLabel(PrintWriter out, MultiNode node) {
        out.print("S\t" + getNodeId(node) + "\t" + node.sequence);
        long coverage = 0;
        for (Map<String, Integer> graph : graphs) {
            for (int i = 0; i + k <= node.sequence.length(); i++) {
                String kmer = node.sequence.substring(i, i + k);
                Integer cov = graph.get(StringUtils.normalizeDna(kmer));
                coverage += cov == null ? 0 : cov;
            }
        }
        out.print("\tLN:i:" + (node.sequence.length()) + "\tKC:i:" + coverage);
        String color = determineColor(node);
        out.println("\tCL:z:" + color + "\tC2:z:" + color);
    }

    private String determineColor(MultiNode node) {
        if (this.graphs.length == 2) {
            if (node.isGeneNode) {
                return COLOR_GREEN;
            }
            switch ((int) node.graphs.size()) {
                case 1:
                    return COLOR_RED;
                case 2:
                    return COLOR_BLUE;
                default:
                    return COLOR_BLACK;
            }
        } if (this.graphs.length == 3) {
            if (node.isGeneNode) {
                return COLOR_GREEN;
            }
            switch ((int) node.mask) {
                case 1:
                    return COLOR_RED;
                case 2:
                    return COLOR_BLUE;
                case 3:
                    return "#ff00ff";
                case 4:
                    return "#ffff00";
                case 5:
                    return "#ffaa00";
                case 6:
                    return "#00ffff";
                default:
                    return COLOR_BLACK;
            }
        } else {
            if (node.isGeneNode) {
                return COLOR_GREEN;
            }
            int value = 256 * node.graphs.size() / graphs.length;
            return String.format("#%02X%02X%02X", value, value, value);
        }
    }

}
