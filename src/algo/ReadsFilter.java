package algo;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.formats.QualityFormat;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

public class ReadsFilter implements Runnable {
    private OneSequenceCalculator calc;
    private File file;
    private final String outputPrefix;
    private final int readsNumber;
    private int k;
    private int percentFiltration;
    private final Logger logger;

    public ReadsFilter(File file, OneSequenceCalculator calc, String outputPrefix, int readsNumber, int k, int percentFiltration, Logger logger) {
        this.file = file;
        this.calc = calc;
        this.outputPrefix = outputPrefix;
        this.readsNumber = readsNumber;
        this.k = k;
        this.percentFiltration = percentFiltration;
        this.logger = logger;
    }

    private void reader() throws ExecutionFailedException {
        //File outputCutReads = new File(outputPrefix + "/cutReads" + readsNumber + ".fastq");
        File outputCutReads = new File(outputPrefix + "/cutReads" + readsNumber + ".fasta");
        outputCutReads.getParentFile().mkdirs();
        NamedSource<Dna> reader = null;
        PrintWriter out = null;
        try {
            out = new PrintWriter(outputCutReads);
            reader = ReadersUtils.readDnaLazy(file);
            ReadersUtils.loadDnaQs(file);
            QualityFormat qualF = ReadersUtils.determineQualityFormat(file);
            NamedSource<DnaQ> quals = ReadersUtils.readDnaQLazy(file);
            int indexCutRead = 0;
            for (DnaQ qual : quals) {
                String read = qual.toString();
                int len = read.length();
                int kmersInRead = len - k + 1;
                int kmersFiltration = kmersInRead * percentFiltration / 100;
                kmersFiltration = Math.max(1, kmersFiltration);
                int numberCoverageKmers = 0;
                for (int i = 0; i < len - k; i++) {
                    if (calc.isContainedInSubgraph(read.substring(i, i + k))) {
                        numberCoverageKmers++;
                    }
                    if (numberCoverageKmers >= kmersFiltration) {
                        indexCutRead++;
                        /*out.println("@" + indexCutRead + "\n" + read + "\n+");
                        for (int j = 0; j < len; j++) {
                            out.print(qualF.getPhredChar(qual.phredAt(j)));
                        }
                        out.println();*/
                        out.println(">" + indexCutRead + "\n" + read);
                        break;
                    }
                }
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Failed to read from file " + file.getPath());
        } finally {
            if (out != null)
                out.close();
        }

    }

    @Override
    public void run() {
        try {
            reader();
        } catch (ExecutionFailedException e) {
            logger.info(e.getMessage());
        }
    }

}

//--k 3 --reads reads.fastq --seq genes.fasta --maxkmers 100 --coverage 0 -o newfile