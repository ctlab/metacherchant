package algo;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.pairs.UniPair;
import utils.HashFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by -- on 18.03.2019.
 */
public abstract class ReadsFinderInGraph {
    final UniPair<LightDnaQ> pair;
    final int k;
    private final BigLong2ShortHashMap graph;
    private final HashFunction hasher;
    final boolean doCorrection;
    private final double z;

    ReadsFinderInGraph(UniPair<LightDnaQ> pair, int k, BigLong2ShortHashMap graph, HashFunction hasher, boolean doCorrection, double z) {
        this.pair = pair;
        this.k = k;
        this.graph = graph;
        this.hasher = hasher;
        this.doCorrection = doCorrection;
        this.z = z;
    }


    boolean findRead(LightDnaQ dnaQ) {
        if (dnaQ.length() < k) {
            return false;
        }
        List<Short> cov = getCoverage(dnaQ);

        double cov_mean = (double) (cov.stream().mapToInt(Short::shortValue).sum() + cov.get(cov.size() - 1) * (k - 1)) / dnaQ.length();
        double width = (double) (cov.stream().mapToInt(i -> i > 0 ? 1 : 0).sum() + (cov.get(cov.size() - 1) > 0 ? 1 : 0) * (k - 1)) / dnaQ.length();
        double theory_width = getTheoryWidth(cov_mean);

        return !(width < 0.9) && delta(cov_mean, width, theory_width, dnaQ.length());
    }

    List<Short> getCoverage(LightDnaQ dnaQ) {
        List<Short> cov = new ArrayList<>();
        if (k > 31) {
            Dna dna = new Dna(dnaQ);
            for (int i = 0; i + k <= dnaQ.length(); i++) {
                long hash = hasher.hash(dna, i, i + k);
                short tmp = graph.getWithZero(hash);
                if (tmp < 0) {
                    throw new RuntimeException("Kmer count < 0");
                }
                cov.add(tmp);
            }
        } else {
            for (ShortKmer kmer : ShortKmer.kmersOf(dnaQ, k)) {
                short tmp = graph.getWithZero(kmer.toLong());
                if (tmp < 0) {
                    throw new RuntimeException("Kmer count < 0");
                }
                cov.add(tmp);
            }
        }
        return cov;
    }

    double getTheoryWidth(double cov) {
        return 1.0 - Math.exp(-cov);
    }

    // normal approximation interval
    boolean delta(double cov, double width, double theory_width, int length) {
        // double z = 1.96; // probit(1 - (1 - confidence_level) / 2) at confidence_level = 0.95
        // z = 1 at confidence_level ~= 0.68
        double std = z * Math.sqrt(Math.exp(-cov) * (1 - Math.exp(-cov)) / length);

        return width == 1 || (width != 0 && -std <= width - theory_width && width - theory_width <= std);
    }

    // Wilson score interval
    boolean delta2(double cov, double width, double t_width, int length) {
        double z = 1.96; // probit(1 - (1 - confidence_level) / 2) at confidence_level = 0.95
        // z = 1 at confidence_level ~= 0.68

        double z2 = Math.pow(z, 2);
        double theory_width = (t_width + z2 / (2 * length)) / (1 + z2 / length);

        double p = Math.exp(-cov);
        double std = z / (1 + z2 / length) * Math.sqrt(p * (1 - p) / length + z2 / (4 * length * length));

        return width == 1 || (width != 0 && -std <= width - theory_width && width - theory_width <= std);
    }

    boolean findReadWithCorrection(LightDnaQ dnaQ) {
        if (dnaQ.length() < k) {
            return false;
        }

        List<Integer> badPos = new ArrayList<>();
        for (int i = 0; i < dnaQ.length(); i++) {
            if (dnaQ.phredAt(i) < 10) {
                badPos.add(i);
            }
        }

        if (badPos.size() > 1) {
            // Do not correct read which has more than 1 error
            //not_found_reads.add(dnaQ);
            return findRead(dnaQ);
        } else if (badPos.size() == 0) {
            return findRead(dnaQ);
        } else {
            // Try to add as found all 4 reads
            boolean found = false;
            for (int nuc = 0; nuc < 4; nuc++) {
                DnaQ correctedDnaQ = new DnaQ(dnaQ);
                correctedDnaQ.setNuc(badPos.get(0), nuc);

                List<Short> cov = getCoverage(correctedDnaQ);

                double cov_mean = (double) (cov.stream().mapToInt(Short::shortValue).sum() + cov.get(cov.size() - 1) * (k - 1)) / correctedDnaQ.length();
                double width = (double) (cov.stream().mapToInt(i -> i > 0 ? 1 : 0).sum() + (cov.get(cov.size() - 1) > 0 ? 1 : 0) * (k - 1)) / correctedDnaQ.length();
                double theory_width = getTheoryWidth(cov_mean);

                if (!(width < 0.9) && delta(cov_mean, width, theory_width, correctedDnaQ.length())) {
                    found = true;
                    break;
                    // TODO: add all four reads or the best (nearest to the theoretical width) ?
                }
            }
            return found;
        }
    }
}
