package tools;

import algo.SingleNode;
import io.IOUtils;
import io.LargeKIOUtils;
import io.writers.GFAWriter;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.tools.io.LazyDnaQReaderTool;
import ru.ifmo.genetics.utils.KmerUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;
import utils.FNV1AHash;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;

import static algo.SingleNode.Color.*;
import static utils.StringUtils.normalizeDna;

/**
 * Created by -- on 20.08.2019.
 */
public class FMTVisualiser extends Tool {
    public static final String NAME = "fmt-visualiser";
    public static final String DESCRIPTION = "Outputs graphs in .gfa format showing the results of FMT classification";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> donorFiles = addParameter(new FileMVParameterBuilder("donor-files")
            .mandatory()
            .withShortOpt("donor")
            .withDescription("file with paired donor metagenomic reads")
            .create());

    public final Parameter<File[]> beforeFiles = addParameter(new FileMVParameterBuilder("before-files")
            .mandatory()
            .withShortOpt("before")
            .withDescription("file with paired pre-FMT recipient metagenomic reads")
            .create());

    public final Parameter<File[]> afterFiles = addParameter(new FileMVParameterBuilder("after-files")
            .mandatory()
            .withShortOpt("after")
            .withDescription("file with paired post-FMT recipient metagenomic reads")
            .create());

    public final Parameter<String> hashFunction = addParameter(new StringParameterBuilder("hash")
            .withDescription("hash function to use: poly or fnv1a")
            .withDefaultValue("poly")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .optional()
            .withShortOpt("o")
            .withDescription("directory to output found reads")
            .withDefaultValue(workDir.append("graph"))
            .create());

    public final Parameter<File> inputDir = addParameter(new FileParameterBuilder("input-dir")
            .mandatory()
            .withShortOpt("i")
            .withDescription("directory containing output of reads_classifier.sh FMT classification script")
            .create());


    private BigLong2ShortHashMap graph, settle, not_settle, stay, gone, from_donor, from_both, from_before, itself;
    private HashFunction hasher;
    private SingleNode[] nodes;
    private String outputPrefix;
    private HashMap<String, Integer> subgraph;
    private int size;

    private long getKmerKey(String s) {
        if (hasher != null) {
            s = normalizeDna(s);
            return hasher.hash(s);
        } else {
            return KmerUtils.getKmerKey(DnaTools.toLong(new Dna(s)), k.get());
        }
    }

    private void loadDonorGraphs(File[] found, File[] not_found) throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.graph = LargeKIOUtils.loadReads(donorFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.settle = LargeKIOUtils.loadReads(found, k.get(), 0, availableProcessors.get(), logger);
            this.not_settle = LargeKIOUtils.loadReads(not_found, k.get(), 0, availableProcessors.get(), logger);
        } else {
            this.graph = IOUtils.loadReads(donorFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.settle = IOUtils.loadReads(found, k.get(), 0, availableProcessors.get(), logger);
            this.not_settle = IOUtils.loadReads(not_found, k.get(), 0, availableProcessors.get(), logger);
        }
        //logger.info("Hashtable size: " + this.graph.size() + " kmers");
    }

    private void loadBeforeGraphs(File[] found, File[] not_found) throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.graph = LargeKIOUtils.loadReads(beforeFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.stay = LargeKIOUtils.loadReads(found, k.get(), 0, availableProcessors.get(), logger);
            this.gone = LargeKIOUtils.loadReads(not_found, k.get(), 0, availableProcessors.get(), logger);
        } else {
            this.graph = IOUtils.loadReads(beforeFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.stay = IOUtils.loadReads(found, k.get(), 0, availableProcessors.get(), logger);
            this.gone = IOUtils.loadReads(not_found, k.get(), 0, availableProcessors.get(), logger);
        }
        //logger.info("Hashtable size: " + this.graph.size() + " kmers");
    }

    private void loadAfterGraphs(File[] from_donor, File[] from_before, File[] from_both, File[] itself) throws ExecutionFailedException {
        if (k.get() > 31) {
            logger.info("Reading hashes of k-mers instead");
            this.hasher = LargeKIOUtils.hash = determineHashFunction();
            this.graph = LargeKIOUtils.loadReads(afterFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.from_donor = LargeKIOUtils.loadReads(from_donor, k.get(), 0, availableProcessors.get(), logger);
            this.from_before = LargeKIOUtils.loadReads(from_before, k.get(), 0, availableProcessors.get(), logger);
            this.from_both = LargeKIOUtils.loadReads(from_both, k.get(), 0, availableProcessors.get(), logger);
            this.itself = LargeKIOUtils.loadReads(itself, k.get(), 0, availableProcessors.get(), logger);
        } else {
            this.graph = IOUtils.loadReads(afterFiles.get(), k.get(), 0, availableProcessors.get(), logger);
            this.from_donor = IOUtils.loadReads(from_donor, k.get(), 0, availableProcessors.get(), logger);
            this.from_before = IOUtils.loadReads(from_before, k.get(), 0, availableProcessors.get(), logger);
            this.from_both = IOUtils.loadReads(from_both, k.get(), 0, availableProcessors.get(), logger);
            this.itself = IOUtils.loadReads(itself, k.get(), 0, availableProcessors.get(), logger);
        }
        //logger.info("Hashtable size: " + this.graph.size() + " kmers");
    }


    private HashFunction determineHashFunction() {
        if (k.get() <= 31) {
            return null;
        }
        String name = hashFunction.get().toLowerCase();
        if (name.equals("fnv1a")) {
            logger.info("Using FNV1a hash function");
            return new FNV1AHash();
        } else {
            logger.info("Using default polynomial hash function");
            return new PolynomialHash();
        }
    }

    @Override
    protected void cleanImpl() {
        graph = null;
        hasher = null;
        settle = null;
        not_settle = null;
        stay = null;
        gone = null;
        from_donor = null;
        from_both = null;
        from_before = null;
        itself = null;
        nodes = null;
        subgraph = null;
    }

    private String toStr(long key) {
        StringBuilder r = new StringBuilder();
        for (int i = 0; i < k.get(); ++i) {
            r.insert(0, DnaTools.NUCLEOTIDES[(int) key & 3]);
            key = key >>> 2;
        }
        return normalizeDna(r.toString());
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        outputPrefix = outputDir.get().getPath();
        String inputPrefix = inputDir.get().getPath() + "/";
        // Drawing donor graph
        {

            File[] settle_files = {new File(inputPrefix + "settle_1.fastq"),
                    new File(inputPrefix + "settle_2.fastq"), new File(inputPrefix + "settle_s.fastq")};
            File[] not_settle_files = {new File(inputPrefix + "not_settle_1.fastq"),
                    new File(inputPrefix + "not_settle_2.fastq"), new File(inputPrefix + "not_settle_s.fastq")};
            loadDonorGraphs(settle_files, not_settle_files);

            logger.info("Loading donor k-mers ...");
            loadKmers(donorFiles.get());

            logger.info("Creating donor image ...");
            createPicture((seq) ->
                    settle.contains(getKmerKey(seq)) && ! not_settle.contains(getKmerKey(seq)) ? GREEN :
                    ! settle.contains(getKmerKey(seq)) && not_settle.contains(getKmerKey(seq)) ? BLUE :
                    settle.contains(getKmerKey(seq)) && not_settle.contains(getKmerKey(seq)) ? GREY :
                    BLACK, "donor");
            cleanImpl();
        }

        // Drawing before graph
        {

            File[] stay_files = {new File(inputPrefix + "stay_1.fastq"),
                    new File(inputPrefix + "stay_2.fastq"), new File(inputPrefix + "stay_s.fastq")};
            File[] gone_files = {new File(inputPrefix + "gone_1.fastq"),
                    new File(inputPrefix + "gone_2.fastq"), new File(inputPrefix + "gone_s.fastq")};
            loadBeforeGraphs(stay_files, gone_files);

            logger.info("Loading before k-mers ...");
            loadKmers(beforeFiles.get());

            logger.info("Creating before image ...");
            createPicture((seq) ->
                    stay.contains(getKmerKey(seq)) && ! gone.contains(getKmerKey(seq)) ? GREEN :
                    ! stay.contains(getKmerKey(seq)) && gone.contains(getKmerKey(seq)) ? BLUE :
                    stay.contains(getKmerKey(seq)) && gone.contains(getKmerKey(seq)) ? GREY :
                    BLACK, "before");
            cleanImpl();
        }

        // Drawing after graph
        {

            File[] from_donor_files = {new File(inputPrefix + "came_from_donor_1.fastq"),
                    new File(inputPrefix + "came_from_donor_2.fastq"),
                    new File(inputPrefix + "came_from_donor_s.fastq")};
            File[] from_before_files = {new File(inputPrefix + "came_from_before_1.fastq"),
                    new File(inputPrefix + "came_from_before_2.fastq"),
                    new File(inputPrefix + "came_from_before_s.fastq")};
            File[] from_both_files = {new File(inputPrefix + "came_from_both_1.fastq"),
                    new File(inputPrefix + "came_from_both_2.fastq"),
                    new File(inputPrefix + "came_from_both_s.fastq")};
            File[] itself_files = {new File(inputPrefix + "came_itself_1.fastq"),
                    new File(inputPrefix + "came_itself_2.fastq"),
                    new File(inputPrefix + "came_itself_s.fastq")};
            loadAfterGraphs(from_donor_files, from_before_files, from_both_files, itself_files);

            logger.info("Loading after k-mers ...");
            loadKmers(afterFiles.get());

            logger.info("Creating after image ...");
            createPicture((seq) ->
                    from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? RED :
                    ! from_donor.contains(getKmerKey(seq)) && from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLUE :
                    ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                        && from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? GREEN :
                    ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && itself.contains(getKmerKey(seq)) ? YELLOW :
                    ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLACK :
                    GREY, "after");
            cleanImpl();
        }
    }

    private void loadKmers(File[] files) throws ExecutionFailedException {
        subgraph = new HashMap<>();
        if (k.get() > 31) {
            LazyDnaQReaderTool dnaQReader = new LazyDnaQReaderTool();
            for (File file : files) {
                dnaQReader.fileIn.set(file);
                dnaQReader.simpleRun();
                NamedSource<? extends DnaQ> source = dnaQReader.dnaQsSourceOut.get();
                for (DnaQ dnaQ : source) {
                    for (int i = 0; i + k.get() <= dnaQ.length(); i++) {
                        String key = dnaQ.substring(i, i + k.get()).toString();
                        subgraph.put(normalizeDna(key), (int) graph.get(getKmerKey(key)));
                    }
                }
            }
        } else {
            graph.entryIterator().forEachRemaining((v) -> subgraph.put(toStr(v.getKey()), (int) v.getValue()));
        }
    }

    private void createPicture(Function<String, SingleNode.Color> getNodeColor, String name) {
        logger.debug("Initializing structures for creating picture ...");
        initializeStructures(getNodeColor);
        logger.debug("Merging vertices for creating picture ...");
        doMerge();
        logger.debug("Outputing merged vertices ...");
        outputNodeSequences(outputPrefix, nodes, name);
        logger.debug("Drawing image ...");
        {
            GFAWriter writer = new GFAWriter(k.get(), outputPrefix, nodes, subgraph, name);
            writer.print();
        }

    }

    private void initializeStructures(Function<String, SingleNode.Color> getNodeColor) {
        Map<String, List<SingleNode>> nodeByKmer = new HashMap<>();
        this.size = subgraph.size() * 2;
        logger.debug("Number of k-mers: " + this.size);
        nodes = new SingleNode[size];
        {
            int id = 0;
            Iterator<Map.Entry<String, Integer>> iter = subgraph.entrySet().iterator();
            while (iter.hasNext()) {
                String seq = iter.next().getKey();
                String rc = DnaTools.reverseComplement(seq);
                SingleNode.Color color = getNodeColor.apply(seq);
                nodes[id] = new SingleNode(seq, id, color);
                nodes[id + 1] = new SingleNode(rc, id + 1, color);
                nodes[id].rc = nodes[id + 1];
                nodes[id + 1].rc = nodes[id];

                id += 2;
            }
        }
        logger.debug("Nodes are created");
        for (int i = 0; i < size; i++) {
            String key = nodes[i].sequence.substring(0, k.get() - 1);
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
                    if (other.neighbors.size() != 1 || nodes[i].color != other.color) {
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


    private void outputNodeSequences(String outputPrefix, SingleNode[] nodes, String name) {
        try {
            File output = new File(outputPrefix + "/" + name + "_seqs.fasta");
            output.getParentFile().mkdirs();
            PrintWriter out = new PrintWriter(output);
            for (int i = 0; i < nodes.length; i++) {
                if (!nodes[i].deleted && nodes[i].id < nodes[i].rc.id && nodes[i].sequence.length() >= 1) {
                    out.print("> ");
                    out.print("Id" + (nodes[i].id + 1) + " ");
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

    private void checkLabels(String a, String b) {
        if (!a.substring(a.length() - (k.get() - 1)).equals(b.substring(0, k.get() - 1))) {
            throw new AssertionError("Labels should be merged, but can not: " + a + " and " + b);
        }
    }

    private String mergeLabels(String a, String b) {
        checkLabels(a, b);
        return a + b.substring(k.get() - 1);
    }

    public FMTVisualiser() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new FMTVisualiser().mainImpl(args);
    }
}
