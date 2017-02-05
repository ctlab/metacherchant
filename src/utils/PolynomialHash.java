package utils;

import ru.ifmo.genetics.dna.Dna;

public class PolynomialHash extends AbstractHashFunction {
    @Override
    public long hash(String s) {
        long hashFW = 1, hashRC = 1;
        for (int i = 0; i < s.length(); i++) {
            hashFW *= 5;
            hashRC *= 5;
            hashFW += nucId[s.charAt(i)];
            hashRC += 3 ^ nucId[s.charAt(s.length() - i - 1)];
        }
        return Math.min(hashFW, hashRC);
    }

    @Override
    public long hash(Dna dna, int start, int end) {
        long hashFW = 1, hashRC = 1;
        for (int i = start; i < end; i++) {
            hashFW *= 5;
            hashRC *= 5;
            hashFW += dna.nucAt(i);
            hashRC += 3 ^ dna.nucAt(end - 1 - (i - start));
        }
        return Math.min(hashFW, hashRC);
    }
}
