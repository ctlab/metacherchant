package algo;

import java.util.ArrayList;
import java.util.List;

public class MultiNode {
    public String sequence;
    public int id;
    public boolean isGeneNode;
    public boolean deleted;
    public MultiNode rc;
    public List<MultiNode> neighbors;
    public long mask;

    public MultiNode(String sequence, int id, boolean isGeneNode) {
        this.sequence = sequence;
        this.id = id;
        this.isGeneNode = isGeneNode;

        this.deleted = false;
        this.neighbors = new ArrayList<MultiNode>();
    }

    public void addFile(int pos) {
        mask |= 1L << pos;
    }
}
