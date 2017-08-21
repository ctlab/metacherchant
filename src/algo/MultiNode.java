package algo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class MultiNode {
    public String sequence;
    public int id;
    public boolean isGeneNode;
    public boolean deleted;
    public MultiNode rc;
    public List<MultiNode> neighbors;
    public Set<Integer> graphs;

    public MultiNode(String sequence, int id, boolean isGeneNode) {
        this.sequence = sequence;
        this.id = id;
        this.isGeneNode = isGeneNode;

        this.deleted = false;
        this.neighbors = new ArrayList<MultiNode>();
        this.graphs = new HashSet<Integer>();
    }

    public void addGraph(int pos) {
        graphs.add(pos);
    }
}
