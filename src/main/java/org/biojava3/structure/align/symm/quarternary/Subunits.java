
package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point3d;

/**
 *
 * @author Peter
 */
public class Subunits {
    private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
    private List<Point3d[]> cbCoords = new ArrayList<Point3d[]>();
    private List<Point3d> originalCenters = new ArrayList<Point3d>(0);
    private List<Point3d> centers = new ArrayList<Point3d>(0);
    private List<Integer> sequenceClusterIds = new ArrayList<Integer>(0);
    private Point3d centroid;
    MomentsOfInertia momentsOfInertia = new MomentsOfInertia();

    public Subunits(List<Point3d[]> caCoords, List<Point3d[]> cbCoords, List<Integer> sequenceClusterIds) {
        this.caCoords = caCoords;
        this.cbCoords = cbCoords;
        this.sequenceClusterIds = sequenceClusterIds;
    }

    public List<Point3d[]> getTraces() {
        return caCoords;
    }
    
    public List<Point3d[]> getCBCoords() {
        return cbCoords;
    }

    public int getSubunitCount() {
        run();
        return centers.size();
    }
    
    public List<Integer> getSequenceClusterIds() {
    	return sequenceClusterIds;
    }
    
    public int getCalphaCount() {
    	int count = 0;
    	for (Point3d[] trace: caCoords) {
    		count += trace.length;
    	}
    	return count;
    }

    public List<Point3d> getCenters() {
        run();
        return centers;
    }

    public List<Point3d> getOrignalCenters() {
        run();
        return originalCenters;
    }

    public Point3d getCentroid() {
        run();
        return centroid;
    }

    public MomentsOfInertia getMomentsOfInertia() {
    	run();
    	return momentsOfInertia;
    }
    private void run() {
        if (centers.size() > 0) {
            return;
        }
        calcOriginalCenters();
        calcCentroid();
        calcCenters();
        calcMomentsOfIntertia();
    }

    private void calcOriginalCenters() {
        for (Point3d[] trace: caCoords) {
            Point3d com = SuperPosition.centroid(trace);
            originalCenters.add(com);
        }
    }

    private void calcCentroid() {
        Point3d[] orig = originalCenters.toArray(new Point3d[originalCenters.size()]);
        centroid = SuperPosition.centroid(orig);
    }

    private void calcCenters() {
        for (Point3d p: originalCenters) {
            Point3d c = new Point3d(p);
            c.sub(centroid);
            centers.add(c);
        }
    }
    
    public Point3d getLowerBound() {
    	Point3d lower = new Point3d();
    	for (Point3d p: centers) {
    		if (p.x < lower.x) {
    			lower.x = p.x;
    		}
    		if (p.y < lower.y) {
    			lower.y = p.y;
    		}
    		if (p.z < lower.z) {
    			lower.z = p.z;
    		}
    	}
    	return lower;
    }
    
    public Point3d getUpperBound() {
    	Point3d upper = new Point3d();
    	for (Point3d p: centers) {
    		if (p.x > upper.x) {
    			upper.x = p.x;
    		}
    		if (p.y > upper.y) {
    			upper.y = p.y;
    		}
    		if (p.z > upper.z) {
    			upper.z = p.z;
    		}
    	}
    	return upper;
    }
    
    private void calcMomentsOfIntertia() {
    	for (Point3d[] trace: caCoords) {
    		for (Point3d p: trace) {
    			momentsOfInertia.addPoint(p, 1.0f);
    		}
    	}
    }
}
