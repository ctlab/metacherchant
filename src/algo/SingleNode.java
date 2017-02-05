package algo;

import java.util.ArrayList;
import java.util.List;

public class SingleNode {
    public String sequence;
    public int id;
    public boolean isGeneNode;
    public boolean deleted;
    public SingleNode rc;
    public List<SingleNode> neighbors;


    public SingleNode(String sequence, int id, boolean isGeneNode) {
        this.sequence = sequence;
        this.id = id;
        this.isGeneNode = isGeneNode;

        this.deleted = false;
        this.neighbors = new ArrayList<SingleNode>();
    }
}
