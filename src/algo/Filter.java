package algo;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.formats.QualityFormat;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.*;
import java.util.Scanner;

/**
 * Created by -- on 07.03.2018.
 */
public class Filter implements Runnable{
    private final String outputPrefix;
    private final int readsNumber;
    private final Logger logger;
    private final Integer maxThreads;

    public Filter(String outputPrefix, int readsNumber, Logger logger, Integer maxThreads) {
        this.outputPrefix = outputPrefix;
        this.readsNumber = readsNumber;
        this.logger = logger;
        this.maxThreads = maxThreads;
    }

    private void runFilter() throws ExecutionFailedException {
        try {
            String[] command = {"blastn", "-db", outputPrefix + "/dbReads",
                    "-task", "blastn-short",
                    "-query", outputPrefix + "/" + readsNumber + ".fasta",
                    "-out", outputPrefix + "/" + readsNumber + ".out",
                    "-num_threads", maxThreads.toString(),
                    "-outfmt", "6 qaccver length pident"};
            ProcessBuilder procBuilder = new ProcessBuilder(command);

            procBuilder.redirectErrorStream(true);
            Process process = procBuilder.start();
            InputStream stdout = process.getInputStream();
            InputStreamReader isrStdout = new InputStreamReader(stdout);
            BufferedReader brStdout = new BufferedReader(isrStdout);
            String line;
            while ((line = brStdout.readLine()) != null) {
                logger.info(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void run() {
        try {
            runFilter();
        } catch (ExecutionFailedException e) {
            logger.info(e.getMessage());
        }
    }
}
