package algo;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.readers.FastqReader;
import ru.ifmo.genetics.io.formats.QualityFormat;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Tool;

public class ReadsFilter implements Runnable {
    private OneSequenceCalculator calc;
    private File[] files;
    private int k;
    private int procentFiltration;

    public ReadsFilter(File[] files, OneSequenceCalculator calc, int k, int procentFiltration) {
        this.files = files;
        this.calc = calc;
        this.k = k;
        this.procentFiltration = procentFiltration;
    }

    public void reader() throws ExecutionFailedException {
        if (files.length == 1) {
            File file = files[0];
            File outputCutReads = new File(calc.outputPrefix + "/cutReads.fastq");
            outputCutReads.getParentFile().mkdirs();
            NamedSource<Dna> reader = null;
            try {
                PrintWriter out;
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
                    int kmersFiltration = kmersInRead * procentFiltration / 100;
                    kmersFiltration = Math.max(1, kmersFiltration);
                    int numberCoverageKmers = 0; 
                    for (int i = 0; i < len - k; i++) {
                        if (calc.isContainedInSubgraph(read.substring(i, i + k))) {
                            numberCoverageKmers++;
                        }
                        if (numberCoverageKmers >= kmersFiltration) {
                            indexCutRead ++;
                            out.println("@" + indexCutRead + "\n" + read + "\n+");
                            for (int j = 0; j < len; j++) {
                                out.print(qualF.getPhredChar(qual.phredAt(j)));
                            }
                            out.println();
                            break;
                        }
                    }
                }
                out.close();
            } catch (IOException e) {
                throw new ExecutionFailedException("Failed to read from file " + file.getPath());
            }
        } 
    }

    @Override
    public void run() {
        try {
            reader();
        } catch (ExecutionFailedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        // TODO Auto-generated method stub

    }

}

//--k 3 --reads reads.fastq --seq genes.fasta --maxkmers 100 --coverage 0 -o newfile