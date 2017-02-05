package utils;

import static ru.ifmo.genetics.dna.DnaTools.NUCLEOTIDES;
import static ru.ifmo.genetics.dna.DnaTools.reverseComplement;

public class StringUtils {

    public static String[] leftNeighbors(String kmer) {
        String[] neighbors = new String[4];
        for (int i = 0; i < 4; i++) {
            neighbors[i] = NUCLEOTIDES[i] + (kmer.substring(0, kmer.length() - 1));
        }
        return neighbors;
    }

    public static String[] rightNeighbors(String kmer) {
        String[] neighbors = new String[4];
        for (int i = 0; i < 4; i++) {
            neighbors[i] = kmer.substring(1) + NUCLEOTIDES[i];
        }
        return neighbors;
    }

    public static String[] allNeighbors(String kmer) {
        String[] neighbors = new String[8];
        String[] left = leftNeighbors(kmer), right = rightNeighbors(kmer);
        for (int i = 0; i < 4; i++) {
            neighbors[2 * i] = left[i];
            neighbors[2 * i + 1] = right[i];
        }
        return neighbors;
    }

    public static String normalizeDna(String s) {
        String rc = reverseComplement(s);
        if (s.compareTo(rc) < 0) {
            return s;
        } else {
            return rc;
        }
    }

    public static String shortenLabel(String label, int k) {
        if (label.length() >= 2 * k) {
            return label.substring(0, k) + "..." + label.substring(label.length() - k) + " (length=" + label.length() + ")";
        } else {
            return label;
        }
    }
}
