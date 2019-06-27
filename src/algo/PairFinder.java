package algo;

import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.pairs.UniPair;
import utils.HashFunction;

import java.util.Queue;

/**
 * Created by -- on 11.03.2019.
 */
public class PairFinder extends ReadsFinderInGraph implements Runnable {
    private final Queue<UniPair<LightDnaQ>> both_found;
    private final Queue<UniPair<LightDnaQ>> first_found;
    private final Queue<UniPair<LightDnaQ>> second_found;
    private final Queue<UniPair<LightDnaQ>> both_not_found;

    public PairFinder(UniPair<LightDnaQ> pair, int k, BigLong2ShortHashMap graph, HashFunction hasher,
                      Queue<UniPair<LightDnaQ>> both_found, Queue<UniPair<LightDnaQ>> first_found,
                      Queue<UniPair<LightDnaQ>> second_found, Queue<UniPair<LightDnaQ>> both_not_found, Boolean doCorrection,
                      double z, double found_threshold) {
        super(pair, k, graph, hasher, doCorrection, z, found_threshold);
        this.both_found = both_found;
        this.first_found = first_found;
        this.second_found = second_found;
        this.both_not_found = both_not_found;
    }

    @Override
    public void run() {
        boolean found_1;
        boolean found_2;
        if (doCorrection) {
            found_1 = findReadWithCorrection(pair.first);
            found_2 = findReadWithCorrection(pair.second);
        } else {
            found_1 = findRead(pair.first);
            found_2 = findRead(pair.second);
        }

        if (pair.second.length() == 0) {
            found_2 = !found_1;
        }

        if (found_1 && found_2) {
            both_found.add(pair);
        } else if (found_1) {
            first_found.add(pair);
        } else if (found_2) {
            second_found.add(pair);
        } else {
            both_not_found.add(pair);
        }
    }
}
