package io.graph;

import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class DeBruijnGraphUtils {
    public static Map<String,Integer> loadGraph(File input) throws ExecutionFailedException {
        BufferedReader br;
        Map<String, Integer> graph = new HashMap<String, Integer>();
        try {
            br = new BufferedReader(new FileReader(input));
            while (br.ready()) {
                String[] tokens = br.readLine().split(" ");
                int depth = Integer.parseInt(tokens[1]);
                graph.put(tokens[0], depth);
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load graph from file " + input.getPath());
        }
        return graph;
    }
}
