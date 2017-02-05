package utils;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;

public interface HashFunction {
    long hash(String s);
    long hash(Dna dna, int start, int end);
}
