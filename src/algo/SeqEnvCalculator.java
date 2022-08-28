package algo;

import io.writers.GFAWriter;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.KmerUtils;
import utils.HashFunction;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;

import static io.writers.GFAWriter.GENE_LABEL_SUFFIX;
import static utils.StringUtils.allNeighbors;
import static utils.StringUtils.normalizeDna;
import static utils.StringUtils.shortenLabel;

/**
 * Created by -- on 22.01.2020.
 */
public class SeqEnvCalculator implements Runnable {
    private final String sequence;
    private final int k;
    private final String outputPrefix;
    private final HashFunction hasher;
    private final BigLong2ShortHashMap graph;
    private final Logger logger;
    private final Function<String, SingleNode.Color> getNodeColor;
    private final String name;
    private final TerminationMode termMode;


    private final Map<String, Integer> subgraph;

    private int size;
    private SingleNode[] nodes;
    private boolean fail = false;

    public SeqEnvCalculator(String sequence, int k, String outputPrefix, HashFunction hasher, BigLong2ShortHashMap graph,
                            Logger logger, Function<String, SingleNode.Color> getNodeColor, String name, TerminationMode termMode) {
        this.sequence = sequence;
        this.k = k;
        this.outputPrefix = outputPrefix;
        this.hasher = hasher;
        this.graph = graph;
        this.logger = logger;
        this.getNodeColor = getNodeColor;
        this.name = name;
        this.termMode = termMode;

        this.subgraph = new HashMap<>();
    }

    @Override
    public void run() {
        logger.info("Finding environment for sequence " + shortenLabel(sequence, k));
        runBfs();
        if (fail) {
            logger.info("Could not find any k-mers of the target gene in the input, halting.");
            return;
        }
        extendEnvironment();

        createPicture();
    }

    private void runBfs() {
        List<String> queue = new ArrayList<String>();
        Map<String, Integer> distanceToKmer = new HashMap<String, Integer>();

        for (int i = 0; i + k <= sequence.length(); i++) {
            String kmer = sequence.substring(i, i + k);
            if (graph.getWithZero(getKmerKey(kmer)) > 0) {
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
            String[] neighbors = allNeighbors(kmer);
            for (String neighbor : neighbors) {
                if (graph.getWithZero(getKmerKey(neighbor)) > 0) {
                    if (termMode.allowsAddition(distanceToKmer, neighbor, distance + 1)) {
                        queue.add(neighbor);
                        distanceToKmer.put(neighbor, distance + 1);
                    }
                }
            }
        }

        for (String kmer : distanceToKmer.keySet()) {
            addToSubgraph(kmer);
        }
    }

    private long getKmerKey(String s) {
        if (hasher != null) {
            s = normalizeDna(s);
            return hasher.hash(s);
        } else {
            return KmerUtils.getKmerKey(DnaTools.toLong(new Dna(s)), k);
        }
    }

    private void addToSubgraph(String kmer) {
        subgraph.put(normalizeDna(kmer), (int) graph.getWithZero(getKmerKey(kmer)));
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
                    if (!subgraph.containsKey(normalizeDna(neighbor)) && graph.getWithZero(getKmerKey(neighbor)) > 0) {
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

    private void createPicture() {
        logger.debug("Initializing structures for creating picture ...");
        initializeStructures();
        logger.debug("Merging vertices for creating picture ...");
        doMerge();
        logger.debug("Outputing merged vertices ...");
        outputNodeSequences(outputPrefix, nodes, name);
        logger.debug("Drawing image ...");
        {
            GFAWriter writer = new GFAWriter(k, outputPrefix, nodes, subgraph, name);
            writer.print();
        }
    }

    private void initializeStructures() {
        Map<String, List<SingleNode>> nodeByKmer = new HashMap<>();
        this.size = subgraph.size() * 2;
        nodes = new SingleNode[size];
        {
            int id = 0;
            Iterator<Map.Entry<String, Integer>> iter = subgraph.entrySet().iterator();
            Set<String> geneSeq = new HashSet<>();
            for (int i = 0; i + k <= sequence.length(); i++) {
                String kmer = sequence.substring(i, i + k);
                geneSeq.add(kmer);
            }
            while (iter.hasNext()) {
                String seq = iter.next().getKey();
                String rc = DnaTools.reverseComplement(seq);
                SingleNode.Color color = getNodeColor.apply(seq);
                boolean isGeneNode = geneSeq.contains(seq) || geneSeq.contains(rc);
                nodes[id] = new SingleNode(seq, id, color, isGeneNode);
                nodes[id + 1] = new SingleNode(rc, id + 1, color, isGeneNode);
                nodes[id].rc = nodes[id + 1];
                nodes[id + 1].rc = nodes[id];

                id += 2;
            }
        }
        logger.debug("Nodes are created");
        for (int i = 0; i < size; i++) {
            String key = nodes[i].sequence.substring(0, k - 1);
            if (!nodeByKmer.containsKey(key)) {
                nodeByKmer.put(key, new ArrayList<>());
            }
            nodeByKmer.get(key).add(nodes[i]);
        }
        logger.debug("Nodes neighbors initialized");
        for (int i = 0; i < size; i++) {
            String lastK = nodes[i].sequence.substring(1);
            if (nodeByKmer.containsKey(lastK)) {
                nodes[i].rc.neighbors.addAll(nodeByKmer.get(lastK));
            }
        }
        logger.debug("Nodes neighbors added");
    }

    private void doMerge() {
        while (true) {
            boolean acted = false;
            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].neighbors.size() == 1) {
                    SingleNode other = nodes[i].neighbors.get(0);
                    if (other.neighbors.size() != 1 || nodes[i].color != other.color || nodes[i].isGeneNode != other.isGeneNode) {
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

    private void mergeNodes(SingleNode firstPlus, SingleNode secondMinus) {
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

    private void checkLabels(String a, String b) {
        if (!a.substring(a.length() - (k - 1)).equals(b.substring(0, k - 1))) {
            throw new AssertionError("Labels should be merged, but can not: " + a + " and " + b);
        }
    }

    private String mergeLabels(String a, String b) {
        checkLabels(a, b);
        return a + b.substring(k - 1);
    }

    private void outputNodeSequences(String outputPrefix, SingleNode[] nodes, String name) {
        try {
            File output = new File(outputPrefix + "/" + name + "_seqs.fasta");
            output.getParentFile().mkdirs();
            PrintWriter out = new PrintWriter(output);
            for (int i = 0; i < nodes.length; i++) {
                if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id && nodes[i].sequence.length() >= 1) {
                    out.print("> ");
                    out.print("Id" + getNodeId(nodes[i]) + " ");
                    out.print("Length:" + nodes[i].sequence.length() + " ");
                    out.print("Neighbors:" + getNeighborIds(nodes[i]));
                    out.println();
                    out.println(nodes[i].sequence);
                }
            }
            out.close();
        } catch (IOException e) {
            logger.info(e.getMessage());
        }
    }

    private String getNodeId(SingleNode node) {
        return (Math.min(node.id, node.rc.id) + 1) + (node.isGeneNode ? GENE_LABEL_SUFFIX : "");
    }

    private Set<Integer> getNeighborIds(SingleNode SingleNode) {
        Set<Integer> result = new TreeSet<>();
        for (SingleNode neighbour : SingleNode.neighbors) {
            result.add(Math.min(neighbour.id, neighbour.rc.id) + 1);
        }
        for (SingleNode neighbour : SingleNode.rc.neighbors) {
            result.add(Math.min(neighbour.id, neighbour.rc.id) + 1);
        }
        result.remove(Math.min(SingleNode.id, SingleNode.rc.id) + 1);
        return result;
    }
}
