package utils;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;

public class FNV1AHash extends AbstractHashFunction {

    private static long FNV_OFFSET_BASIS = -3750763034362895579L; // 14695981039346656037L;
    private static long FNV_PRIME = 1099511628211L;

    private static final int DEFAULT_READ_LENGTH = 101;

    private long[] bufFW, bufRC;

    public FNV1AHash() {
        this.bufFW = new long[DEFAULT_READ_LENGTH];
        this.bufRC = new long[DEFAULT_READ_LENGTH];
    }

    @Override
    public long hash(String s) {
        long hashFW = FNV_OFFSET_BASIS, hashRC = FNV_OFFSET_BASIS;
        for (int i = 0; i < s.length(); i++) {
            hashFW ^= nucId[s.charAt(i)];
            hashRC ^= 3 ^ nucId[s.charAt(s.length() - i - 1)];
            hashFW *= FNV_PRIME;
            hashRC *= FNV_PRIME;
        }
        return Math.min(hashFW, hashRC);
    }

    @Override
    public long hash(Dna dna, int start, int end) {
        long hashFW = FNV_OFFSET_BASIS, hashRC = FNV_OFFSET_BASIS;
        for (int i = start; i < end; i++) {
            hashFW ^= dna.nucAt(i);
            hashRC ^= 3 ^ dna.nucAt(start + end - 1 - i);
            hashFW *= FNV_PRIME;
            hashRC *= FNV_PRIME;
        }
        return Math.min(hashFW, hashRC);
    }
}

