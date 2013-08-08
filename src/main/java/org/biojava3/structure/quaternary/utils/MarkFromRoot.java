package org.biojava3.structure.quaternary.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 *
 * @author Peter
 */
public class MarkFromRoot<V> {
    
    private Graph<V> graph;
    private V root = null;
    private V start = null;
    private boolean[] encounter = null;
    private List<V> markedVertices = new ArrayList<V>();
    
    /** Creates a new instance of ComponentFinder */
    public MarkFromRoot() {
    }
    
    public void setGraph(Graph<V> graph) {
        this.graph = graph;
        encounter = new boolean[graph.size()];
    }
    
    public void setRootVertex(V root) {
       this.root = root;
    }
    
    public void setStartVertex(V start) {
        this.start = start;
    }
    
    public List<V> getMarkedVertices() {
        mark();
        return markedVertices;
    }
    
    private void mark() {
        // initially all vertices are marked as not encountered
        markedVertices.clear();
        Arrays.fill(encounter, false);
        
        // mark root atom so we can't traverse 
        // in the direction of the root atom
        encounter[graph.indexOf(root)] = true;

        // mark start atom and add to marked vertices list
        int iStart = graph.indexOf(start);
        encounter[iStart] = true;
        markedVertices.add(start);
        
        // traverse graph from start atom
        traverse(iStart, markedVertices);
    }
    
    private void traverse(int start, List<V> fragment) {
        // mark all neighbors and add to fragment (recursively)
        List<Integer> neighbors = graph.getNeighborIndices(start);
        
        // if there are no neighbors, return
        if (neighbors.size() == 0) return;
        
        for (int neighbor: neighbors) {
            if (! encounter[neighbor]) {
                fragment.add(graph.getVertex(neighbor));
                encounter[neighbor] = true;
                traverse(neighbor, fragment);
            }
        }
    }
    
}
