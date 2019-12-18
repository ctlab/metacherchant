package tools;

import algo.KmerEnvCalculator;
import algo.SingleNode;
import io.IOUtils;
import io.LargeKIOUtils;
import io.LargeKmerLoader;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import static algo.SingleNode.Color.*;
import static algo.SingleNode.Color.GREY;
import static utils.StringUtils.normalizeDna;

/**
 * Created by -- on 18.12.2019.
 */
public class FMTVisualizer extends Tool {
    public static final String NAME = "fmt-visualizer";
    public static final String DESCRIPTION = "Outputs disconnected components of graphs in .gfa format showing the results of FMT classification";

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

    public final Parameter<String> extension = addParameter(new StringParameterBuilder("ext")
            .mandatory()
            .withShortOpt("ext")
            .withDescription("extension of output files of reads_classifier.sh FMT classification script")
            .create());


    private BigLong2ShortHashMap graph, settle, not_settle, stay, gone, from_donor, from_both, from_before, itself;
    private HashFunction hasher;
    private List<DnaQ> sequences;
    private String outputPrefix;

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
        try {
            this.sequences = ReadersUtils.loadDnaQs(donorFiles.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Could not load sequences from " + donorFiles.get()[0].getPath());
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
        try {
            this.sequences = ReadersUtils.loadDnaQs(beforeFiles.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Could not load sequences from " + donorFiles.get()[0].getPath());
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
        try {
            this.sequences = ReadersUtils.loadDnaQs(afterFiles.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Could not load sequences from " + donorFiles.get()[0].getPath());
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
        sequences = null;
        stay = null;
        gone = null;
        from_donor = null;
        from_both = null;
        from_before = null;
        itself = null;
    }

    private Function<String, SingleNode.Color> getDonorColorNode() {
        return (seq) -> settle.contains(getKmerKey(seq)) && ! not_settle.contains(getKmerKey(seq)) ? GREEN :
                        ! settle.contains(getKmerKey(seq)) && not_settle.contains(getKmerKey(seq)) ? BLUE :
                                settle.contains(getKmerKey(seq)) && not_settle.contains(getKmerKey(seq)) ? GREY :
                                        BLACK;
    }

    private Function<String, SingleNode.Color> getBeforeColorNode() {
        return (seq) -> stay.contains(getKmerKey(seq)) && ! gone.contains(getKmerKey(seq)) ? GREEN :
                        ! stay.contains(getKmerKey(seq)) && gone.contains(getKmerKey(seq)) ? BLUE :
                                stay.contains(getKmerKey(seq)) && gone.contains(getKmerKey(seq)) ? GREY :
                                        BLACK;
    }

    private Function<String, SingleNode.Color> getAfterColorNode() {
        return (seq) -> from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? RED :
                        ! from_donor.contains(getKmerKey(seq)) && from_before.contains(getKmerKey(seq))
                                && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLUE :
                                ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                        && from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? GREEN :
                                        ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                                && ! from_both.contains(getKmerKey(seq)) && itself.contains(getKmerKey(seq)) ? YELLOW :
                                                ! from_donor.contains(getKmerKey(seq)) && ! from_before.contains(getKmerKey(seq))
                                                        && ! from_both.contains(getKmerKey(seq)) && ! itself.contains(getKmerKey(seq)) ? BLACK :
                                                        GREY;
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        outputPrefix = outputDir.get().getPath();
        String inputPrefix = inputDir.get().getPath() + "/";
        // Drawing donor graph
        {

            logger.info("Loading donor reads ...");
            File[] settle_files = {new File(inputPrefix + "settle_1." + extension.get()),
                    new File(inputPrefix + "settle_2." + extension.get()),
                    new File(inputPrefix + "settle_s." + extension.get())};
            File[] not_settle_files = {new File(inputPrefix + "not_settle_1." + extension.get()),
                    new File(inputPrefix + "not_settle_2." + extension.get()),
                    new File(inputPrefix + "not_settle_s." + extension.get())};
            loadDonorGraphs(settle_files, not_settle_files);

            logger.info("Creating donor image ...");
            long comp = 0;
            for (DnaQ dseq : sequences) {
                String seq = dseq.toString();
                for (int i = 0; i + k.get() <= seq.length(); i++) {
                    String kmer = seq.substring(i, i + k.get());
                    if (graph.get(getKmerKey(kmer)) > 0) {
                        new KmerEnvCalculator(kmer, k.get(), outputPrefix + "/donor", this.hasher, graph, logger,
                                getDonorColorNode(), "comp" + comp).run();
                        comp++;
                    }

                }
            }
            cleanImpl();
        }

        // Drawing before graph
        {
            logger.info("Loading before reads ...");
            File[] stay_files = {new File(inputPrefix + "stay_1." + extension.get()),
                    new File(inputPrefix + "stay_2." + extension.get()),
                    new File(inputPrefix + "stay_s." + extension.get())};
            File[] gone_files = {new File(inputPrefix + "gone_1." + extension.get()),
                    new File(inputPrefix + "gone_2." + extension.get()),
                    new File(inputPrefix + "gone_s." + extension.get())};
            loadBeforeGraphs(stay_files, gone_files);

            logger.info("Creating before image ...");
            long comp = 0;
            for (DnaQ dseq : sequences) {
                String seq = dseq.toString();
                for (int i = 0; i + k.get() <= seq.length(); i++) {
                    String kmer = seq.substring(i, i + k.get());
                    if (graph.get(getKmerKey(kmer)) > 0) {
                        new KmerEnvCalculator(kmer, k.get(), outputPrefix + "/before", this.hasher, graph, logger,
                                getBeforeColorNode(), "comp" + comp).run();
                        comp++;
                    }

                }
            }
            cleanImpl();
        }

        // Drawing after graph
        {
            logger.info("Loading after reads ...");
            File[] from_donor_files = {new File(inputPrefix + "came_from_donor_1." + extension.get()),
                    new File(inputPrefix + "came_from_donor_2." + extension.get()),
                    new File(inputPrefix + "came_from_donor_s." + extension.get())};
            File[] from_before_files = {new File(inputPrefix + "came_from_before_1." + extension.get()),
                    new File(inputPrefix + "came_from_before_2." + extension.get()),
                    new File(inputPrefix + "came_from_before_s." + extension.get())};
            File[] from_both_files = {new File(inputPrefix + "came_from_both_1." + extension.get()),
                    new File(inputPrefix + "came_from_both_2." + extension.get()),
                    new File(inputPrefix + "came_from_both_s." + extension.get())};
            File[] itself_files = {new File(inputPrefix + "came_itself_1." + extension.get()),
                    new File(inputPrefix + "came_itself_2." + extension.get()),
                    new File(inputPrefix + "came_itself_s." + extension.get())};
            loadAfterGraphs(from_donor_files, from_before_files, from_both_files, itself_files);

            logger.info("Creating after image ...");
            long comp = 0;
            for (DnaQ dseq : sequences) {
                String seq = dseq.toString();
                for (int i = 0; i + k.get() <= seq.length(); i++) {
                    String kmer = seq.substring(i, i + k.get());
                    if (graph.get(getKmerKey(kmer)) > 0) {
                        new KmerEnvCalculator(kmer, k.get(), outputPrefix + "/after", this.hasher, graph, logger,
                                getAfterColorNode(), "comp" + comp).run();
                        comp++;
                    }

                }
            }
            cleanImpl();
        }
    }

    public FMTVisualizer() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new FMTVisualizer().mainImpl(args);
    }
}
