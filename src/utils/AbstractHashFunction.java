package utils;

import ru.ifmo.genetics.dna.DnaTools;

public abstract class AbstractHashFunction implements HashFunction {

    static byte[] nucId = new byte[256];
    static {
        for (byte i = 0; i < 4; i++) {
            nucId[DnaTools.NUCLEOTIDES[i]] = i;
        }
    }


}
