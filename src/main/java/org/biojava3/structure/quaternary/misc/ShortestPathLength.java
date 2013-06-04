package org.biojava3.structure.quaternary.misc;
import java.util.List;

import org.biojava3.structure.utils.SimpleGraph;

/**
 *
 * @author Peter
 */
public class ShortestPathLength<V> {
    
    private V temp;
    private List<Integer> stack;
    private int pathCounter;
    private SimpleGraph<V> graph;
    private V startVertex;
    private int startIndex;
    private int[] pathLength;
    private int[] currentElements;
    private int[] nextElements;
    private int currentElementCount;
    private int nextElementCount;
    private int radius;
    /**
     * Creates a new instance of ShortestPathLength
     */
    
    public ShortestPathLength(SimpleGraph<V> graph) {
        this.graph = graph;
        this.startVertex = null;
        pathLength = new int[graph.size()];
        currentElements = new int[graph.size()];
        currentElementCount = 0;
        nextElements = new int[graph.size()];
        nextElementCount = 0;
        for (int i = 0; i < pathLength.length; i++) pathLength[i] = -1;
        int pathCounter = 0;
        radius = graph.size()-1;
    }
    
    public void setRadius(int radius) {
        this.radius = radius;
    }
    
    public void setStartVertex(V startVertex) {
        this.startVertex = startVertex;
        startIndex = graph.indexOf(startVertex);
        search();
    }
    
    public int getPathLength(V vertex) {
        int index = graph.indexOf(vertex);
        if (index == -1) return index;
        return pathLength[index];
    }
    
    public int[][] getDistanceMatrix() {
        int [][] triangle = new int[graph.size()][];
        for (int i = 0, n = graph.size(); i < n; i++) {
            triangle[i] = new int[i + 1];
            startIndex = i;
            search();
            for (int j = 0; j < i + 1; j++) {
                triangle[i][j] = pathLength[j];
            }
        }
        return triangle;
    }
    
    private void search() {
        pathCounter = 0;
        currentElements[0] = startIndex;
        currentElementCount = 1;
        
        for (int i = 0; i < pathLength.length; i++) pathLength[i] = -1;
        pathLength[currentElements[0]] = 0;
        
        while (currentElementCount > 0 && pathCounter < radius) {
            pathCounter++;
            nextElementCount = 0;
            for (int i = 0; i < currentElementCount; i++) {
                nextIteration(currentElements[i]);
            }
            int[] tmp = currentElements;
            currentElements = nextElements;
            currentElementCount = nextElementCount;
            nextElements = tmp;
        }
        
    }
    
    private void nextIteration(int root) {
        for (Integer i: graph.getNeighborIndices(root)) {
            if (pathLength[i] == -1) {
                pathLength[i] = pathCounter;
                nextElements[nextElementCount] = i;
                nextElementCount++;
            }
        }
    }
    
}
