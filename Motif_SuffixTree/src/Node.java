/*
Created on April 28, 2017

@author: Hannie
*/

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

class Node {

    private int[] data;
    private int lastIdx = 0;
    private static final int START_SIZE = 0;
    private static final int INCREMENT = 1;
    private final Map<Character, Edge> edges;
    private Node suffix;
    
    //# of different results that are stored in this node and in underlying ones 
    private int resultCount = -1;

    Node() {
        edges = new EdgeBag();
        suffix = null;
        data = new int[START_SIZE];
    }

    //get children
    Collection<Integer> getData() {
        return getData(-1);
    }

    Collection<Integer> getData(int numElements) {
        Set<Integer> ret = new HashSet<Integer>();
        for (int num : data) {
            ret.add(num);
            if (ret.size() == numElements) {
                return ret;
            }
        }
        // need to get more matches from child nodes. This is what may waste time
        for (Edge e : edges.values()) {
            if (-1 == numElements || ret.size() < numElements) {
                for (int num : e.getDest().getData()) {
                    ret.add(num);
                    if (ret.size() == numElements) {
                        return ret;
                    }
                }
            }
        }
        return ret;
    }

    /**
     * Adds the given <tt>index</tt> to the set of indexes associated with <tt>this</tt>
     */
    void addRef(int index) {
        if (contains(index)) {
            return;
        }

        addIndex(index);

        // add this reference to all the suffixes as well
        Node iter = this.suffix;
        while (iter != null) {
            if (iter.contains(index)) {
                break;
            }
            iter.addRef(index);
            iter = iter.suffix;
        }

    }

    private boolean contains(int index) {
        int low = 0;
        int high = lastIdx - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            int midVal = data[mid];

            if (midVal < index)
            low = mid + 1;
            else if (midVal > index)
            high = mid - 1;
            else
            return true;
        }
        return false;	// Crap! Java 1.8 does not work with Utils.BinarySearch
    }

//    Computes the number of results that are stored on this node and on its
//    children, and caches the result.
    protected int computeAndCacheCount() {
        computeAndCacheCountRecursive();
        return resultCount;
    }

    private Set<Integer> computeAndCacheCountRecursive() {
        Set<Integer> ret = new HashSet<Integer>();
        for (int num : data) {
            ret.add(num);
        }
        for (Edge e : edges.values()) {
            for (int num : e.getDest().computeAndCacheCountRecursive()) {
                ret.add(num);
            }
        }

        resultCount = ret.size();
        return ret;
    }


//     Returns the number of results that are stored on this node and on its
//     children.
//     Should be called after having called computeAndCacheCount.

    public int getResultCount() throws IllegalStateException {
        if (-1 == resultCount) {
            throw new IllegalStateException();
        }

        return resultCount;
    }

    void addEdge(char ch, Edge e) {
        edges.put(ch, e);
    }

    Edge getEdge(char ch) {
        return edges.get(ch);
    }

    Map<Character, Edge> getEdges() {
        return edges;
    }

    Node getSuffix() {
        return suffix;
    }

    void setSuffix(Node suffix) {
        this.suffix = suffix;
    }

    private void addIndex(int index) {
        if (lastIdx == data.length) {
            int[] copy = new int[data.length + INCREMENT];
            System.arraycopy(data, 0, copy, 0, data.length);
            data = copy;
        }
        data[lastIdx++] = index;
    }
}
