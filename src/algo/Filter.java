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
    private final String database;
    private final String databasePath;
    private final String outputPrefix;
    private final int readsNumber;
    private final Logger logger;
    private final Integer maxThreads;

    public Filter(String database, String databasePath, String outputPrefix, int readsNumber, Logger logger, Integer maxThreads) {
        this.database = database;
        this.databasePath = databasePath;
        this.outputPrefix = outputPrefix;
        this.readsNumber = readsNumber;
        this.logger = logger;
        this.maxThreads = maxThreads;
    }

    private void runFilter() throws ExecutionFailedException {
        try {
            String[] command = {"blastn", "-db", databasePath + "/" + database,
                    "-query", outputPrefix + "cutReads" + readsNumber + ".fasta",
                    "-out", outputPrefix + "filter" + readsNumber + ".txt",
                    "-num_threads", maxThreads.toString()};
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
        File reads = new File(outputPrefix + "cutReads" + readsNumber + ".fasta");
        File outputFiltReads = new File(outputPrefix + "filtReads" + readsNumber + ".fasta");
        outputFiltReads.getParentFile().mkdirs();
        NamedSource<Dna> reader = null;
        PrintWriter out = null;

        try {
            Scanner filtered = new Scanner(new File(outputPrefix + "filter" + readsNumber + ".txt"));
            out = new PrintWriter(outputFiltReads);
            reader = ReadersUtils.readDnaLazy(reads);
            ReadersUtils.loadDnaQs(reads);
            QualityFormat qualF = ReadersUtils.determineQualityFormat(reads);
            NamedSource<DnaQ> quals = ReadersUtils.readDnaQLazy(reads);
            int indexFiltRead = 1;
            int indexRead = 0;
            for (DnaQ qual : quals) {
                indexRead++;
                String read = qual.toString();
                while (filtered.hasNextLine() && !filtered.nextLine().startsWith("Query= " + indexRead)) {
                }
                filtered.nextLine();
                filtered.nextLine();
                filtered.nextLine();
                String line = filtered.nextLine();
                if (line.startsWith("Sequences producing significant alignments:")) {
                    out.println(">" + (indexFiltRead++) + "\n" + read);
                } else if (filtered.nextLine().startsWith("***** No hits found *****")) {
                    continue;
                } else {
                    throw new ExecutionFailedException("Illegal result of blastn search found in file " +
                            (outputPrefix + "filtReads" + readsNumber + ".txt"));
                }
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Failed to read from file " + reads.getPath());
        } finally {
            if (out != null)
                out.close();
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
