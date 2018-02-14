package algo;

import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by -- on 13.02.2018.
 */
public class AssemblerCalculator implements Runnable {
    private final String assembler;
    private final String assemblerPath;
    private final String outputPrefix;
    private final Logger logger;

    public AssemblerCalculator(String assembler, String assemblerPath, String outputPrefix, Logger logger) {
        this.assembler = assembler;
        this.assemblerPath = assemblerPath;
        this.outputPrefix = outputPrefix;
        this.logger = logger;
    }

    private void runAssembler() {
        if (assembler.equals("spades")) {
            try {
                String[] command = {"python", assemblerPath + "/spades.py", "--12",
                        outputPrefix + "cutReads.fastq", "-o", outputPrefix + "out_spades"};
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
                // or better take assembly_graph.gfa?
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if (assembler.equals("megahit")) {
            try {
                String[] commandA = {assemblerPath + "/megahit", "--12",
                        outputPrefix + "cutReads.fastq", "-o", outputPrefix + "out_megahit"};
                String[] commandR = {"mv", outputPrefix + "out_megahit/" + "final.contigs.fa",
                        outputPrefix + "out_megahit/" + "final.contigs.fasta"};
                ProcessBuilder procBuilder = new ProcessBuilder(commandA);
                ProcessBuilder procBuilderR = new ProcessBuilder(commandR);

                procBuilder.redirectErrorStream(true);
                Process process = procBuilder.start();
                InputStream stdout = process.getInputStream();
                InputStreamReader isrStdout = new InputStreamReader(stdout);
                BufferedReader brStdout = new BufferedReader(isrStdout);

                String line;
                while ((line = brStdout.readLine()) != null) {
                    logger.info(line);
                }

                procBuilderR.redirectErrorStream(true);
                Process processB = procBuilderR.start();
                InputStream stdoutB = processB.getInputStream();
                InputStreamReader isrStdoutB = new InputStreamReader(stdoutB);
                BufferedReader brStdoutB = new BufferedReader(isrStdoutB);
                String lineB;
                while ((lineB = brStdoutB.readLine()) != null) {
                    logger.info(lineB);
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void run() {
        runAssembler();
    }
}
