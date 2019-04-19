package algo;

import ru.ifmo.genetics.dna.LightDnaQ;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.pairs.UniPair;
import tools.TripleReadsClassifier;
import utils.HashFunction;

import java.util.List;
import java.util.Map;
import java.util.Queue;

/**
 * Created by -- on 18.03.2019.
 */
public class TripleFinder2 extends ReadsFinderInGraph implements Runnable {
    private final Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_1;
    private final Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_2;
    private final Queue<UniPair<LightDnaQ>> both_found;
    private final Queue<UniPair<LightDnaQ>> both_half_found;
    private final Queue<UniPair<LightDnaQ>> both_not_found;
    private final Queue<LightDnaQ> s_found;
    private final Queue<LightDnaQ> s_half_found;
    private final Queue<LightDnaQ> s_not_found;

    public TripleFinder2(UniPair<LightDnaQ> pair, int k, BigLong2ShortHashMap graph, HashFunction hasher,
                         Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_1,
                         Map<String, TripleReadsClassifier.FindResult> isFoundInGraphOne_2, boolean doCorrection,
                         Queue<UniPair<LightDnaQ>> both_found, Queue<UniPair<LightDnaQ>> both_half_found,
                         Queue<UniPair<LightDnaQ>> both_not_found, Queue<LightDnaQ> s_found,
                         Queue<LightDnaQ> s_half_found, Queue<LightDnaQ> s_not_found, double z) {
        super(pair, k, graph, hasher, doCorrection, z);
        this.isFoundInGraphOne_1 = isFoundInGraphOne_1;
        this.isFoundInGraphOne_2 = isFoundInGraphOne_2;
        this.both_found = both_found;
        this.both_half_found = both_half_found;
        this.both_not_found = both_not_found;
        this.s_found = s_found;
        this.s_half_found = s_half_found;
        this.s_not_found = s_not_found;
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

        TripleReadsClassifier.FindResult res_1, res_2;
        if (found_1 && TripleReadsClassifier.FindResult.FOUND.equals(isFoundInGraphOne_1.get(pair.first.toString()))) {
            res_1 =  TripleReadsClassifier.FindResult.FOUND;
        } else if (found_1 || TripleReadsClassifier.FindResult.FOUND.equals(isFoundInGraphOne_1.get(pair.first.toString())) ||
                (getWidth(pair.first) >= 0.4 && TripleReadsClassifier.FindResult.HALF_FOUND.equals(isFoundInGraphOne_1.get(pair.first.toString())))) {
            res_1 = TripleReadsClassifier.FindResult.HALF_FOUND;
        } else {
            res_1 = TripleReadsClassifier.FindResult.NOT_FOUND;
        }

        if (found_2 && TripleReadsClassifier.FindResult.FOUND.equals(isFoundInGraphOne_2.get(pair.second.toString()))) {
            res_2 =  TripleReadsClassifier.FindResult.FOUND;
        } else if (found_2 || TripleReadsClassifier.FindResult.FOUND.equals(isFoundInGraphOne_2.get(pair.second.toString())) ||
                (getWidth(pair.second) >= 0.4 && TripleReadsClassifier.FindResult.HALF_FOUND.equals(isFoundInGraphOne_2.get(pair.second.toString())))) {
            res_2 = TripleReadsClassifier.FindResult.HALF_FOUND;
        } else {
            res_2 = TripleReadsClassifier.FindResult.NOT_FOUND;
        }


        if (res_1 == TripleReadsClassifier.FindResult.FOUND && res_2 == TripleReadsClassifier.FindResult.FOUND) {
            both_found.add(pair);
        } else if (res_1 == TripleReadsClassifier.FindResult.FOUND && res_2 == TripleReadsClassifier.FindResult.HALF_FOUND) {
            s_found.add(pair.first);
            s_half_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.FOUND && res_2 == TripleReadsClassifier.FindResult.NOT_FOUND) {
            s_found.add(pair.first);
            s_not_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.HALF_FOUND && res_2 == TripleReadsClassifier.FindResult.FOUND) {
            s_half_found.add(pair.first);
            s_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.HALF_FOUND && res_2 == TripleReadsClassifier.FindResult.HALF_FOUND) {
           both_half_found.add(pair);
        } else if (res_1 == TripleReadsClassifier.FindResult.HALF_FOUND && res_2 == TripleReadsClassifier.FindResult.NOT_FOUND) {
            s_half_found.add(pair.first);
            s_not_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.NOT_FOUND && res_2 == TripleReadsClassifier.FindResult.FOUND) {
            s_not_found.add(pair.first);
            s_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.NOT_FOUND && res_2 == TripleReadsClassifier.FindResult.HALF_FOUND) {
            s_not_found.add(pair.first);
            s_half_found.add(pair.second);
        } else if (res_1 == TripleReadsClassifier.FindResult.NOT_FOUND && res_2 == TripleReadsClassifier.FindResult.NOT_FOUND) {
            both_not_found.add(pair);
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
