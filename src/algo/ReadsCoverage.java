package algo;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by -- on 14.03.2018.
 */
public class ReadsCoverage implements Runnable{
    private final String outputPrefix;
    private final String workPrefix;
    private final int readsNumber;
    private final Logger logger;

    public ReadsCoverage(String outputPrefix, String workPrefix, int readsNumber, Logger logger) {
        this.outputPrefix = outputPrefix;
        this.workPrefix = workPrefix;
        this.readsNumber = readsNumber;
        this.logger = logger;
    }

    private void runFilter() throws ExecutionFailedException {
        try {
            StringBuilder input = new StringBuilder("\"");
            for (int i = 0; i < readsNumber; i++) {
                input.append(outputPrefix).append("cutReads").append(i).append(".fasta ");
            }
            input.append("\"");
            String[] command = {"makeblastdb", "-in", input.toString(),
                    "-parse_seqids", "-dbtype", "nucl",
                    "-out", workPrefix + "db/" + "dbReads"};
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
