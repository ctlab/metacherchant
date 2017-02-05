package algo;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class TerminationMode {
    public enum TerminationModeType {
        MAX_KMERS,
        MAX_RADIUS
    }

    private final List<TerminationModeType> types;
    private final List<Integer> thresholds;

    public TerminationMode() {
        types = new ArrayList<TerminationModeType>();
        thresholds = new ArrayList<Integer>();
    }

    public TerminationMode(TerminationModeType type, int threshold) {
        this();
        addRestriction(type, threshold);
    }

    public void addRestriction(TerminationModeType type, int threshold) {
        types.add(type);
        thresholds.add(threshold);
    }

    public boolean allowsAddition(Map<String, Integer> distanceToKmer, String kmer, int newDistance) {
        if (distanceToKmer.containsKey(kmer)) {
            return false;
        }
        for (int i = 0; i < types.size(); i++) {
            int threshold = thresholds.get(i);
            switch (types.get(i)) {
                case MAX_KMERS:
                    if (distanceToKmer.size() >= threshold) return false;
                    break;
                case MAX_RADIUS:
                    if (newDistance > threshold) return false;
                    break;
            }
        }
        return true;
    }
}
