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
    boolean visited;
    boolean changed;

    public Color color;

    public enum Color {RED, GREEN, BLUE, GREY, YELLOW, BLACK}


    public SingleNode(String sequence, int id, boolean isGeneNode) {
        this.sequence = sequence;
        this.id = id;
        this.isGeneNode = isGeneNode;
        this.color = null;

        this.deleted = false;
        this.visited = false;
        this.changed = false;
        this.neighbors = new ArrayList<SingleNode>();
    }

    public SingleNode(String sequence, int id, Color color) {
        this.sequence = sequence;
        this.id = id;
        this.isGeneNode = false;
        this.color = color;

        this.deleted = false;
        this.visited = false;
        this.changed = false;
        this.neighbors = new ArrayList<SingleNode>();
    }
}
