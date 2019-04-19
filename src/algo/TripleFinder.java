package algo;

import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.pairs.UniPair;
import tools.TripleReadsClassifier;
import utils.HashFunction;

import java.util.List;
import java.util.Map;

/**
 * Created by -- on 18.03.2019.
 */
public class TripleFinder extends ReadsFinderInGraph implements Runnable {
    private final Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_1;
    private final Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_2;

    public TripleFinder(UniPair<LightDnaQ> pair, int k, BigLong2ShortHashMap graph, HashFunction hasher,
                        Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_1,
                        Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_2, boolean doCorrection, double z) {
        super(pair, k, graph, hasher, doCorrection, z);
        this.isFoundInGraphOne_1 = isFoundInGraphOne_1;
        this.isFoundInGraphOne_2 = isFoundInGraphOne_2;
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

        if (found_1) {
            isFoundInGraphOne_1.put(pair.first.toString(), TripleReadsClassifier.FindResult.FOUND);
        } else if (getWidth(pair.first) >= 0.4) {
            isFoundInGraphOne_1.put(pair.first.toString(), TripleReadsClassifier.FindResult.HALF_FOUND);
        } else {
            isFoundInGraphOne_1.put(pair.first.toString(), TripleReadsClassifier.FindResult.NOT_FOUND);
        }
        if (found_2) {
            isFoundInGraphOne_2.put(pair.second.toString(), TripleReadsClassifier.FindResult.FOUND);
        } else if (getWidth(pair.second) >= 0.4) {
            isFoundInGraphOne_2.put(pair.second.toString(), TripleReadsClassifier.FindResult.HALF_FOUND);
        } else {
            isFoundInGraphOne_2.put(pair.second.toString(), TripleReadsClassifier.FindResult.NOT_FOUND);
        }

    }

    private double getWidth(LightDnaQ dnaQ) {
        if (dnaQ.length() < k) {
            return 0;
        }
        List<Short> cov = getCoverage(dnaQ);
        return (double) (cov.stream().mapToInt(i -> i > 0 ? 1 : 0).sum() + (cov.get(cov.size() - 1) > 0 ? 1 : 0) * (k - 1)) / dnaQ.length();
    }
}
