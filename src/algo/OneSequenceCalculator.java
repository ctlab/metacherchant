package algo;

import io.writers.GFAWriter;
import io.writers.TSVWriter;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.KmerUtils;
import utils.HashFunction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import static utils.StringUtils.*;

public class OneSequenceCalculator implements Runnable {
    private final String sequence;
    private final String outputPrefix;
    private final int minOccurences;
    private final int k;
    private final HashFunction hasher;
    private final BigLong2ShortHashMap reads;
    private final boolean bothDirections;
    private final Logger logger;
    private final int chunkLength;
    private final TerminationMode termMode;
    private final boolean trimPaths;

    private final Map<String, Integer> subgraph;

    private int size;
    private SingleNode[] nodes;
    private boolean fail = false;

    public OneSequenceCalculator(String sequence, int k, int minOccurences, String outputPrefix, HashFunction hasher, BigLong2ShortHashMap reads, Logger logger, boolean bothDirections, int chunkLength, TerminationMode termMode, boolean trimPaths) {
        this.sequence = sequence;
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.minOccurences = minOccurences;
        this.hasher = hasher;
        this.reads = reads;
        this.logger = logger;
        this.bothDirections = bothDirections;
        this.chunkLength = chunkLength;
        this.termMode = termMode;
        this.trimPaths = trimPaths;

        this.subgraph = new HashMap<String, Integer>();
    }

    private long getKmerKey(String s) {
        if (hasher != null) {
            s = normalizeDna(s);
            return hasher.hash(s);
        } else {
            return KmerUtils.getKmerKey(DnaTools.toLong(new Dna(s)), k);
        }
    }

    @Override
    public void run() {
        logger.info("Finding environment for sequence " + shortenLabel(sequence, k));
        buildEnvironment();
        if (fail) {
            logger.info("Could not find any k-mers of the target gene in the input, halting.");
            return;
        }
        extendEnvironment();
//        printGene();
        printEnvironment();
        createPicture();
    }

    private void printGene() {
        File output = new File(outputPrefix + "/gene.fasta");
        PrintWriter out;
        try {
            out = new PrintWriter(output);
            out.println(sequence);
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void buildEnvironment() {
        if (bothDirections) {
            runBfs(0);
        } else {
            runBfs(-1);
            runBfs(1);
        }
    }

    private void addToSubgraph(String kmer) {
        subgraph.put(normalizeDna(kmer), (int) reads.get(getKmerKey(kmer)));
    }

    boolean isContainedInSubgraph(String kmer) {
        return subgraph.containsKey(normalizeDna(kmer));
    }

    void runBfs(int dir) { // -1 - backward, +1 - forward, 0 - both
        List<String> queue = new ArrayList<String>();
        Map<String, Integer> distanceToKmer = new HashMap<String, Integer>();
        Set<String> lastKmers = new HashSet<String>();

        for (int i = 0; i + k <= sequence.length(); i++) {
            String kmer = sequence.substring(i, i + k);
            int occs = reads.get(getKmerKey(kmer));
            if (occs >= minOccurences) {
                queue.add(kmer);
                distanceToKmer.put(kmer, 0);
            }
        }
        if (queue.size() == 0) {
            fail = true;
            return;
        }
        int head = 0;
        while (head < queue.size()) {
            String kmer = queue.get(head++);
            int distance = distanceToKmer.get(kmer);
            String[] neighbors = getNeighborsByDir(dir, kmer);
            for (String neighbor : neighbors) {
                int occs = reads.get(getKmerKey(neighbor));
                if (occs >= minOccurences) {
                    if (termMode.allowsAddition(distanceToKmer, neighbor, distance + 1)) {
                        queue.add(neighbor);
                        distanceToKmer.put(neighbor, distance + 1);
                    } else {
                        lastKmers.add(kmer);
                    }
                }
            }
        }
        if (trimPaths) {
            runTrimPaths(lastKmers, distanceToKmer, dir);
        }
        for (String kmer : distanceToKmer.keySet()) {
            addToSubgraph(kmer);
        }
    }

    private String[] getNeighborsByDir(int dir, String kmer) {
        String[] neighbors = null;
        switch (dir) {
            case -1: {
                neighbors = leftNeighbors(kmer);
                break;
            }
            case 1: {
                neighbors = rightNeighbors(kmer);
                break;
            }
            case 0: {
                neighbors = allNeighbors(kmer);
                break;
            }
        }
        return neighbors;
    }

    private void runTrimPaths(Set<String> lastKmers, Map<String, Integer> distanceToKmer, int dir) {
        List<String> queue = new ArrayList<String>();
        Set<String> visitedKmers = new HashSet<String>();

        for (String kmer : lastKmers) {
            queue.add(kmer);
            visitedKmers.add(kmer);
        }

        int head = 0;
        while (head < queue.size()) {
            String kmer = queue.get(head++);
            String[] neighbors = getNeighborsByDir(-dir, kmer);
            for (String neighbor : neighbors) {
                if (distanceToKmer.containsKey(neighbor) && !visitedKmers.contains(neighbor)) {
                    queue.add(neighbor);
                    visitedKmers.add(neighbor);
                }
            }
        }
        distanceToKmer.keySet().retainAll(visitedKmers);
    }


    private void extendEnvironment() {
        Iterator<Map.Entry<String, Integer>> iter = subgraph.entrySet().iterator();
        Set<String> additions = new HashSet<String>();
        while (iter.hasNext()) {
            String kmer = iter.next().getKey();
            while (true) {
                String[] neighbors = allNeighbors(kmer);
                String cont = null;
                for (String neighbor : neighbors) {
                    if (!isContainedInSubgraph(neighbor) && reads.get(getKmerKey(neighbor)) >= minOccurences) {
                        if (cont == null) {
                            cont = kmer;
                        } else {
                            cont = "";
                        }
                    }
                }
                if (cont != null && !cont.equals("") && !additions.contains(cont)) {
                    additions.add(cont);
                    kmer = cont;
                } else {
                    break;
                }
            }
        }

        logger.info("Extending endings by " + additions.size() + " kmers");
        for (String kmer : additions) {
            addToSubgraph(kmer);
        }
    }

    private void printEnvironment() {
        File outputEnv = new File(outputPrefix + "/graph.txt");
        outputEnv.getParentFile().mkdirs();
        try {
            PrintWriter out;
            out = new PrintWriter(outputEnv);
            for (Map.Entry<String, Integer> entry : subgraph.entrySet()) {
                out.println(entry.getKey() + " " + entry.getValue());
            }
            out.close();
        } catch (FileNotFoundException e) {
            logger.info("Could not write environment to " + outputEnv.getPath());
        }
    }

    void mergeNodes(SingleNode firstPlus, SingleNode secondMinus) {
        // first k-1 symbols of firstPlus coincide with complement of first k-1 symbols of secondMinus
        SingleNode firstMinus = firstPlus.rc, secondPlus = secondMinus.rc;
        String newSeq = mergeLabels(secondPlus.sequence, firstPlus.sequence);
        String newSeqRC = mergeLabels(firstMinus.sequence, secondMinus.sequence);

        secondPlus.sequence = newSeq;
        firstMinus.sequence = newSeqRC;
        secondPlus.rc = firstMinus;
        firstMinus.rc = secondPlus;

        firstPlus.deleted = secondMinus.deleted = true;
    }

    private void createPicture() {
        initializeStructures();
        doMerge();
        outputNodeSequences();

        {
            GFAWriter writer = new GFAWriter(k, outputPrefix, nodes, subgraph);
            writer.print();
        }
        {
            TSVWriter writer = new TSVWriter(k, outputPrefix, nodes, subgraph);
            writer.print();
        }
    }

    void outputNodeSequences() {
        try {
            PrintWriter out = new PrintWriter(outputPrefix + "/seqs.fasta");
            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id && nodes[i].sequence.length() >= chunkLength) {
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
        }
    }

    private Set<Integer> getNeighborIds(SingleNode node) {
        Set<Integer> result = new TreeSet<Integer>();
        for (SingleNode neighbor : node.neighbors) {
            result.add(Math.min(neighbor.id, neighbor.rc.id) + 1);
        }
        for (SingleNode neighbor : node.rc.neighbors) {
            result.add(Math.min(neighbor.id, neighbor.rc.id) + 1);
        }
        result.remove(Math.min(node.id, node.rc.id) + 1);
        return result;
    }

    private void initializeStructures() {
        Map<String, List<SingleNode>> nodeByKmer = new HashMap<String, List<SingleNode>>();
        this.size = subgraph.size() * 2;
        nodes = new SingleNode[size];
        {
            int id = 0;
            Iterator<Map.Entry<String, Integer>> iter = subgraph.entrySet().iterator();
            while (iter.hasNext()) {
                String seq = iter.next().getKey();
                String rc = DnaTools.reverseComplement(seq);
                boolean isGeneNode = sequence.contains(seq) || sequence.contains(rc);
                nodes[id] = new SingleNode(seq, id, isGeneNode);
                nodes[id + 1] = new SingleNode(rc, id + 1, isGeneNode);
                nodes[id].rc = nodes[id + 1];
                nodes[id + 1].rc = nodes[id];

                id += 2;
            }
        }
        for (int i = 0; i < size; i++) {
            String key = nodes[i].sequence.substring(0, k - 1);
            if (!nodeByKmer.containsKey(key)) {
                nodeByKmer.put(key, new ArrayList<SingleNode>());
            }
            nodeByKmer.get(key).add(nodes[i]);
        }
        for (int i = 0; i < size; i++) {
            String lastK = nodes[i].sequence.substring(1);
            if (nodeByKmer.containsKey(lastK)) {
                nodes[i].rc.neighbors.addAll(nodeByKmer.get(lastK));
            }
        }
    }

    private void doMerge() {
        while (true) {
            boolean acted = false;
            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].neighbors.size() == 1) {
                    SingleNode other = nodes[i].neighbors.get(0);
                    if (other.neighbors.size() != 1 || nodes[i].isGeneNode != other.isGeneNode) {
                        continue;
                    }
                    mergeNodes(nodes[i], other);
                    acted = true;
                }
            }
            if (!acted) {
                break;
            }
        }
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
