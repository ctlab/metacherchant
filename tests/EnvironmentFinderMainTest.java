import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import tools.EnvironmentFinderMain;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;


public class EnvironmentFinderMainTest {
    @Rule
    public TemporaryFolder folder = new TemporaryFolder();

    @Test
    public void testHiCSeqPositive() throws Exception {
        String outputPath = folder.getRoot().getPath() + "output";
        String workDirPath = folder.getRoot().getPath() + "workDir";
        EnvironmentFinderMain main = new EnvironmentFinderMain();
        String[] args = {"--k", "31", "--coverage", "5",
                "--reads", "./Hi-C_pipline/example/wgs_reads/wgs_reads_R1.fastq",
                "./Hi-C_pipline/example/wgs_reads/wgs_reads_R2.fastq",
                "--seq", "./Hi-C_pipline/example/seq.fasta",
                "--output", outputPath,
                "--work-dir", workDirPath,
                "--maxradius", "100000", "--bothdirs", "False", "--chunklength", "10", "--merge", "true",
                "--hicseq", "./Hi-C_pipline/example_work_dir/1/selected_reads.fasta"};
        main.mainImpl(args);

        assert (main.hicSeqsFile != null);
        Field hicSequences = EnvironmentFinderMain.class.getDeclaredField("hicSequences");
        hicSequences.setAccessible(true);
        assertEquals (1047, ((List<DnaQ>) hicSequences.get(main)).size());
        File graph = new File(outputPath + "/merged/graph.gfa");
        assert (graph.exists());
        assert (checkGraph(outputPath + "/merged/graph.gfa", 19, 22));
    }

    @Test
    public void testSeq_WrongFileName() throws Exception {
        EnvFinderTest main = new EnvFinderTest();
        String[] args = {"--k", "31", "--coverage", "5",
                "--reads", "./Hi-C_pipline/example/wgs_reads/wgs_reads_R1.fastq",
                "./Hi-C_pipline/example/wgs_reads/wgs_reads_R2.fastq",
                "--seq", "./Hi-C_pipline/example/wrong_name.fasta",
                "--output", "./Hi-C_pipline/example/junite_test_hic/output",
                "--work-dir", "./Hi-C_pipline/example/junite_test_hic/workDir",
                "--maxradius", "100000", "--bothdirs", "False", "--chunklength", "10", "--merge", "true"};

        main.addParam(args);
        ExecutionFailedException expectedException = null;
        try {
            main.loadInput();
        } catch (ExecutionFailedException ex) {
            expectedException = ex;
        } finally {
            //deleteDir(new File("./Hi-C_pipline/example/junite_test_hic/output"));
            //deleteDir(new File("./Hi-C_pipline/example/junite_test_hic/workDir"));
            assert(expectedException != null);
            assertEquals(expectedException.getMessage(), "Could not load sequences from .\\Hi-C_pipline\\example\\wrong_name.fasta");
        }
    }

    @Test
    public void testHiCSeq_WrongFileName() throws Exception {
        EnvFinderTest main = new EnvFinderTest();
        String[] args = {"--k", "31", "--coverage", "5",
                "--reads", "./Hi-C_pipline/example/wgs_reads/wgs_reads_R1.fastq",
                "./Hi-C_pipline/example/wgs_reads/wgs_reads_R2.fastq",
                "--seq", "./Hi-C_pipline/example/seq.fasta",
                "--output", "./Hi-C_pipline/example/junite_test_hic/output",
                "--work-dir", "./Hi-C_pipline/example/junite_test_hic/workDir",
                "--maxradius", "100000", "--bothdirs", "False", "--chunklength", "10", "--merge", "true",
                "--hicseq", "./Hi-C_pipline/example_work_dir/1/wrong_name.fasta"};

        main.addParam(args);
        ExecutionFailedException expectedException = null;
        try {
            main.loadInput();
        } catch (ExecutionFailedException ex) {
            expectedException = ex;
        } finally {
            assert(expectedException != null);
            assertEquals(expectedException.getMessage(), "Could not load Hi-C sequences from .\\Hi-C_pipline\\example_work_dir\\1\\wrong_name.fasta");
        }
    }

    private boolean checkGraph(String graphGFAPath, int expectedNodesCount, int expectedEdgesCount) throws Exception {
        BufferedReader reader;
        int actualNodesCount = 0;
        int actualEdgesCount = 0;
        reader = new BufferedReader(new FileReader(graphGFAPath));
        String line = reader.readLine();
        while (line != null) {
            if (line.split("\t")[0] .equals("S")) {
                actualNodesCount += 1;
            } else if (line.split("\t")[0] .equals("L")) {
                actualEdgesCount += 1;
            }
            line = reader.readLine();
        }
        reader.close();
        actualEdgesCount /= 2;
        return (expectedNodesCount == actualNodesCount && expectedEdgesCount == actualEdgesCount);
    }

    private class EnvFinderTest extends EnvironmentFinderMain {
        public void addParam(String[] args) {

            List<Parameter> allParameters = new ArrayList<>();
            allParameters.addAll(this.inputParameters);
            allParameters.addAll(launchOptions);
            parseArgs(allParameters, args);
        }
    }
}


