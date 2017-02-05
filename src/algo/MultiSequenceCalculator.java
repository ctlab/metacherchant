package algo;

import io.writers.GFAWriterMulti;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.DnaTools;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import static utils.StringUtils.*;

public class MultiSequenceCalculator implements Runnable {
    private final String sequence;
    private final String outputPrefix;
    private final int k;
    private final Logger logger;
    private final Map<String, Integer>[] graphs;

    private int size;
    private MultiNode[] nodes;

    public MultiSequenceCalculator(String sequence, int k, String outputPrefix, Logger logger, Map<String, Integer>[] graphs) {
        this.sequence = sequence;
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.logger = logger;
        this.graphs = graphs;
    }

    @Override
    public void run() {
        logger.info("Combining environments for sequence " + shortenLabel(sequence, k));
        createPicture();
    }

    private void createPicture() {
        initializeStructures();
        doMerge();
        outputNodeSequences();

        {
            GFAWriterMulti writer = new GFAWriterMulti(k, outputPrefix, nodes, graphs);
            writer.print();
        }
    }

    private void initializeStructures() {
        Map<String, MultiNode> nodeByKmer = new HashMap<String, MultiNode>();
        for (Map<String, Integer> graph : graphs) {
            for (String kmer : graph.keySet()) {
                nodeByKmer.put(kmer, null);
                nodeByKmer.put(DnaTools.reverseComplement(kmer), null);
            }
        }
        this.size = nodeByKmer.size();
        nodes = new MultiNode[size];
        {
            int ptr = 0;
            for (String kmer : nodeByKmer.keySet()) {
                String rc = DnaTools.reverseComplement(kmer);
                if (kmer.compareTo(rc) > 0) {
                    continue;
                }

                boolean isGeneNode = sequence.contains(kmer) || sequence.contains(rc);
                nodes[ptr] = new MultiNode(kmer, ptr, isGeneNode);
                nodes[ptr + 1] = new MultiNode(rc, ptr + 1, isGeneNode);
                nodes[ptr].rc = nodes[ptr + 1];
                nodes[ptr + 1].rc = nodes[ptr];

                nodeByKmer.put(nodes[ptr].sequence, nodes[ptr]);
                nodeByKmer.put(nodes[ptr + 1].sequence, nodes[ptr + 1]);

                ptr += 2;
            }
        }
        for (int i = 0; i < graphs.length; i++) {
            Map<String, Integer> subgraph = graphs[i];
            for (String kmer : subgraph.keySet()) {
                MultiNode node = nodeByKmer.get(kmer);
                node.addFile(i);
                node.rc.addFile(i);
            }
        }

        for (int i = 0; i < nodes.length; i++) {
            String seq = nodes[i].sequence;
            for (int j = 0; j < 4; j++) {
                String next = seq.substring(1) + DnaTools.NUCLEOTIDES[j];
                MultiNode neighbor = nodeByKmer.get(next);
                if (neighbor != null) {
                    nodes[i].rc.neighbors.add(neighbor);
                }
            }
        }
    }

    private void doMerge() {
        while (true) {
            boolean acted = false;
            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].neighbors.size() == 1) {
                    MultiNode other = nodes[i].neighbors.get(0);
                    if (canBeMerged(nodes[i], other)) {
                        mergeNodes(nodes[i], other);
                        acted = true;
                    }
                }
            }
            if (!acted) {
                break;
            }
        }
    }

    private boolean canBeMerged(MultiNode first, MultiNode second) {
        return first.neighbors.size() == 1 && second.neighbors.size() == 1 && first.isGeneNode == second.isGeneNode && first.mask == second.mask;
    }

    private void mergeNodes(MultiNode firstPlus, MultiNode secondMinus) {
        // first k-1 symbols of firstPlus coincide with complement of first k-1 symbols of secondMinus
        MultiNode firstMinus = firstPlus.rc, secondPlus = secondMinus.rc;
        String newSeq = mergeLabels(secondPlus.sequence, firstPlus.sequence);
        String newSeqRC = mergeLabels(firstMinus.sequence, secondMinus.sequence);

        secondPlus.sequence = newSeq;
        firstMinus.sequence = newSeqRC;
        secondPlus.rc = firstMinus;
        firstMinus.rc = secondPlus;

        firstPlus.deleted = secondMinus.deleted = true;
    }

    private void outputNodeSequences() {
        try {
            File outputFile = new File(outputPrefix + "/seqs.fasta");
            outputFile.getParentFile().mkdirs();
            PrintWriter out = new PrintWriter(outputFile);

            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id) {
                    out.print("> ");
                    out.print("Id" + (i + 1) + " ");
                    out.print("Length:" + nodes[i].sequence.length() + " ");
                    out.print("Neighbors:" + getNeighborIds(nodes[i]));
                    out.println();
                    out.println(nodes[i].sequence);
                }
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Set<Integer> getNeighborIds(MultiNode node) {
        Set<Integer> result = new TreeSet<Integer>();
        for (MultiNode neighbor : node.neighbors) {
            result.add(Math.min(neighbor.id, neighbor.rc.id) + 1);
        }
        for (MultiNode neighbor : node.rc.neighbors) {
            result.add(Math.min(neighbor.id, neighbor.rc.id) + 1);
        }
        result.remove(Math.min(node.id, node.rc.id) + 1);
        return result;
    }

    private void checkLabels(String a, String b) {
        if (!a.substring(a.length() - (k - 1)).equals(b.substring(0, k - 1))) {
            throw new AssertionError("Labels should be merged, but can not: " + a + " and " + b);
        }
    }

    private String mergeLabels(String a, String b) {
        checkLabels(a, b);
        return a + b.substring(k - 1);
    }
}
