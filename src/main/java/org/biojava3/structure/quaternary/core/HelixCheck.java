package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Color4f;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.geometry.SuperPosition;
import org.biojava3.structure.quaternary.misc.ColorBrewer;
import org.biojava3.structure.quaternary.utils.ComponentFinder;
import org.biojava3.structure.quaternary.utils.Graph;
import org.biojava3.structure.quaternary.utils.MarkFromRoot;
import org.biojava3.structure.utils.SimpleGraph;

// 3J06, 1PFI (part helix, then expanding radius) helical
public class HelixCheck {
	private static final Vector3d Y_AXIS = new Vector3d(0,1,0);
	private static final Vector3d Z_AXIS = new Vector3d(0,0,1);
	private static String N_FOLD_AXIS_COLOR = "red";
	private static double AXIS_SCALE_FACTOR = 1.4;
	
	private Subunits subunits = null;
	private RotationGroup rotationGroup = null; 
	private QuatSymmetryParameters quatSymmetryParameters = null;
	
	private boolean helical = false;
	private Vector3d minBoundary = new Vector3d();
	private Vector3d maxBoundary = new Vector3d();
	private double xzRadiusMax = Double.MIN_VALUE;
	
	private List<HelixParameters> helixList = new ArrayList<HelixParameters>();

	public HelixCheck(Subunits subunits, RotationGroup rotationGroup, QuatSymmetryParameters quatSymmetryParameters) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
		this.quatSymmetryParameters = quatSymmetryParameters;
	}
	
	public boolean isHelical() {
		if (! precheck()) {
			System.out.println("HelixCheck: failed precheck");
			return false;
		}
		System.out.println("HelixCheck: passed precheck");
		run();
		if (helixList.isEmpty()) {
			System.out.println("HelixCheck: no helix found");
			return false;
		}
//		this.defaultHelixParameters = getByLargestInteraction();
//		this.defaultHelixParameters = getByLowestStartNumber();
		HelixParameters lowestTwistAngle = getByLowestTwistAngle();
		if (Math.abs(lowestTwistAngle.getRise()) < 0.1) {
			System.out.println("HelixCheck: helix rise < 0.1");
			return false;
		}
		System.out.println("Helix parameters: " + lowestTwistAngle);
		System.out.println(getDefaultOrientation(lowestTwistAngle));
		System.out.println(colorBySymmetry(lowestTwistAngle));

		Matrix4d transformation = alignHelixAxis(lowestTwistAngle);
		System.out.println(drawHelixAxis(lowestTwistAngle, transformation));
		System.out.println(getJmolLatticeLines(lowestTwistAngle, transformation, 0));
	    HelixParameters largest = getByLargestInteraction();
	    if (largest == lowestTwistAngle) {
	    	largest = getByLowestStartNumber();
	    }
	    System.out.println(getJmolLatticeLines(largest, transformation, 100));
	    System.out.println("draw axes* on; draw latt* on");
	    return true;
	}
	
	/**
	 * Returns a Jmol script to set the default orientation for a structure
	 * @return Jmol script
	 */
	public String getDefaultOrientation(HelixParameters helixParameters) {	
		StringBuilder s = new StringBuilder();
		
		// calculate  orientation
		Quat4d q = new Quat4d();
		q.set(alignHelixAxis(helixParameters));
		
		// set orientation
		s.append("moveto 0 quaternion{");
		s.append(jMolFloat(q.x));
		s.append(",");
		s.append(jMolFloat(q.y));
		s.append(",");
		s.append(jMolFloat(q.z));
		s.append(",");
		s.append(jMolFloat(q.w));
		s.append("};");
		return s.toString();
	}
	
	/**
	 * Returns a Jmol script that colors subunits to highlight the symmetry within a structure
	 * @return Jmol script
	 */
	public String colorBySymmetry(HelixParameters helixParameters) {
		List<Integer> modelNumbers = this.subunits.getModelNumbers();
		List<String> chainIds = this.subunits.getChainIds();

		List<List<Integer>> layerLines = helixParameters.getLayerLines();
		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
		
		if (layerLines.size() > 1) {
			Color4f[] colors = getSymmetryColors(layerLines.size()); 
			for (int i = 0; i < layerLines.size(); i++) {
				Color4f c = colors[i];
				List<String> ids = colorMap.get(c);
				if (ids == null) {
					ids = new ArrayList<String>();
					colorMap.put(c,  ids);
				}
				for (int subunit: layerLines.get(i)) {
					String id = chainIds.get(subunit) + "/" + (modelNumbers.get(subunit)+1);
					ids.add(id);
				}
			}
		} else {
			List<Integer> oneStartLine = layerLines.get(0);
			Color4f[] colors = getSymmetryColors(oneStartLine.size()); 
			for (int i = 0; i < oneStartLine.size(); i++) {
				int subunit = oneStartLine.get(i);
				String id = chainIds.get(subunit) + "/" + (modelNumbers.get(subunit)+1);
				List<String> ids = Collections.singletonList(id);
				Color4f c = colors[i];
				colorMap.put(c,  ids);
			}
		}
		return getJmolColorScript(colorMap);
	}
	
	public String drawHelixAxis(HelixParameters hp, Matrix4d transformation) {
		StringBuilder s = new StringBuilder();
		calcBoundaries(transformation);

		s.append("draw axesHelix");
		s.append(0);
		s.append(" ");
		s.append("line");
		Point3d v1 = new Point3d(hp.getAxisVector());
		v1.scale(AXIS_SCALE_FACTOR*(maxBoundary.y-minBoundary.y)*0.5);
		Point3d v2 = new Point3d(v1);
		v2.negate();
		v1.add(subunits.getCentroid());
		v2.add(subunits.getCentroid());
		s.append(getJmolPoint(v1));
		s.append(getJmolPoint(v2));
		s.append("width 0.5 ");
		s.append(" color ");
		s.append(N_FOLD_AXIS_COLOR);
		s.append(" off;");
		return s.toString();
	};
	
	
//	private Point3d getGeometricCenter(Matrix4d transformation) {
//		Point3d centroid = new Point3d(subunits.getCentroid());
//		transformation.transform(centroid);
//		
//	}
	private static String getJmolColorScript(Map<Color4f, List<String>> map) {
		StringBuilder s = new StringBuilder();
		for (Entry<Color4f, List<String>> entry: map.entrySet()) {
			s.append("select ");
			List<String> ids = entry.getValue();
			for (int i = 0; i < ids.size(); i++) {
				s.append("*:");
				s.append(ids.get(i));
				if (i < ids.size() -1 ) {
					s.append(",");
				} else {
					s.append(";");
				}
			}
			s.append("color cartoon");	
			s.append(getJmolColor(entry.getKey()));
			s.append(";");
			s.append("color atom");
			s.append(getJmolColor(entry.getKey()));
			s.append(";");

		}
		return s.toString();
	}
	
	/**
	 * Returns a color palette for helical structures
	 * @param nColors
	 * @return
	 */
	private static Color4f[] getSymmetryColors(int nColors) {
		Color4f[] colors = ColorBrewer.RdYlBu.getColor4fPalette(nColors);
		return colors;	
	}
	
	private void run() {
		for (Entry<Integer,List<Integer>> entry: getPairs().entrySet()) {
			int contacts = (int)Math.floor(entry.getKey()/1000.0);
			alignPair(entry.getValue(), contacts);
		}
		sortLayerLines();
	}
	
	private void alignPair(List<Integer> pair, int interactions) {
		Point3d[][]	hl = tracePair(pair);

		Matrix4d transformation = SuperPosition.superposeWithTranslation(hl[0], hl[1]);
		double rmsd = SuperPosition.rmsd(hl[0], hl[1]);
		
		System.out.println("-----------------------------------------------------------");
		System.out.println("HelixCheck.helixAlignment() rmsd: " + rmsd);

		List<Integer> permutations = getPermutation(transformation);
		System.out.println("Permutation: " + permutations);
		List<List<Integer>> layerLines = getLayerLines(permutations);
		System.out.println("Layer lines: " + layerLines);
        helixAlignment(layerLines, interactions);
	}
	
	private List<Integer> getPermutation(Matrix4d transformation) {
		//	double rmsdThresholdSq = Math.pow(this.quatSymmetryParameters.getRmsdThreshold(), 2);
		double rmsdThresholdSq = 49;
		List<Point3d> centers = subunits.getOriginalCenters();

		List<Integer> permutations = new ArrayList<Integer>(centers.size());
		double[] dSqs = new double[centers.size()];

		for (Point3d center: centers) {
			Point3d tCenter = new Point3d(center);
			transformation.transform(tCenter);
			int permutation = -1;
			double minDistSq = Double.MAX_VALUE;
			for (int i = 0; i < centers.size(); i++) {
				double dSq = tCenter.distanceSquared(centers.get(i));
				if (dSq < minDistSq && dSq <= rmsdThresholdSq) {
					minDistSq = dSq;
					permutation = i; 
					dSqs[i] = dSq;
				}
			}
            permutations.add(permutation);
		}
		
		// keep only assignments for the shortest distance
		for (int i = 0; i < centers.size(); i++) {
			double minDistSq = Double.MAX_VALUE;
		    // find the minimum distance
			for (int j = 0; j < centers.size(); j++) {
				if (permutations.get(i) == j) {
					if (dSqs[j] < minDistSq) {
						minDistSq = dSqs[j];
					}
				}
			}
			if (dSqs[i] > minDistSq) {
//				permutations.set(i, -1);
			}
		}
		return permutations;
	}
	
	private static List<List<Integer>> getLayerLines(List<Integer> permutations) {
		Graph<Integer> graph = new SimpleGraph<Integer>();
		for (int i = 0; i < permutations.size(); i++) {
			if (permutations.get(i) != -1) {
				graph.addVertex(i);
				graph.addVertex(permutations.get(i));
				graph.addEdge(i, permutations.get(i));
			}
//			System.out.println("createGraph.addEdge: " + pairs.get(i) +"," + pairs.get(i+1));
		}
		ComponentFinder<Integer> finder = new ComponentFinder<Integer>();
		finder.setGraph(graph);
//		System.out.println("Components:" + finder.getComponentCount());
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (List<Integer> c: finder.getComponents()) {
//			System.out.println("input path: " +c);
//			System.out.println("components: " + findPath(graph, c));
			List<Integer> path = findPath(graph, c);
			if (! path.isEmpty()) {
				paths.add(path);
			}
		}
		return paths;
	}
	
	private Point3d[][] tracePair(List<Integer> pair) {		
		List<Point3d[]> traces = this.subunits.getTraces();
		List<Point3d> xList = new ArrayList<Point3d>();
		List<Point3d> yList = new ArrayList<Point3d>();

		System.out.println("HelixCheck.tracePairs: " + pair.get(0) + "-" + pair.get(1));
		for (Point3d p: traces.get(pair.get(0))) {
			xList.add(new Point3d(p));
		}
		for (Point3d p: traces.get(pair.get(1))) {
			yList.add(new Point3d(p));
		}

		Point3d[] x = xList.toArray(new Point3d[xList.size()]);
		Point3d[] y = yList.toArray(new Point3d[yList.size()]);
		Point3d[][] hl = {x, y};
		return hl;
	}
	
	private Map<Integer, List<Integer>> getPairs() {
		List<Point3d[]> traces = this.subunits.getTraces();
		List<Point3d> centers = this.subunits.getOriginalCenters();
		List<Integer> seqClusterIds = this.subunits.getSequenceClusterIds();
		TreeMap<Integer, List<Integer>> layers = new TreeMap<Integer, List<Integer>>();
		for (int i = 0; i < traces.size()-1; i++) {
			int seqClusterI = seqClusterIds.get(i);
			for (int j = i+1; j < traces.size(); j++) {
				if (seqClusterIds.get(j) != seqClusterI) {
					continue;
				}
                Integer contacts = calcContactNumber(traces.get(i), traces.get(j));        
  //              System.out.println("d: " + centers.get(i).distance(centers.get(j)) + " c: " + contacts);
                contacts = (int)Math.round(contacts*0.2)*5;
                if (contacts == 0) {
                	continue;
                }
                double d = centers.get(i).distance(centers.get(j));
                d = (int)Math.round(d*10)/10;
               
                Integer key = contacts*1000 + (int)d;
               
				List<Integer> layer = layers.get(key);
				if (layer == null) {
					System.out.println(i + "," + j + "  d: " + d + " cbin: " + contacts + " key: " + key);
					layer = new ArrayList<Integer>();
					layers.put(key, layer);
					layer.add(i);
					layer.add(j);
				}
			}
		}
		// only keep the top 4 layers
		System.out.println("Layers before: " + layers.size());
		while (layers.size() > 4) {
			Integer firstKey = layers.firstKey();
			layers.remove(firstKey);
		}
		System.out.println("Layers after: " + layers.size());
		return layers;
	}
	
	private static int calcContactNumber(Point3d[] a, Point3d[] b) {
		int contacts = 0;
		for (Point3d pa : a) {
			for (Point3d pb : b) {
				if (pa.distance(pb) < 10) {
					contacts++;
				}
			}
		}
		return contacts;
	}
	
	private static List<Integer> findPath(Graph<Integer> graph, List<Integer> components) {
		Integer root = null;
		for (Integer vertex: graph.getVertices()) {
			if (graph.getValence(vertex) == 1 && components.contains(vertex)) {
				root = vertex;
				break;
			}
		}
		if (root != null) {
		   MarkFromRoot<Integer> marker = new MarkFromRoot<Integer>();
		   marker.setGraph(graph);
		   marker.setRootVertex(root);
		   marker.setStartVertex(root);
		   List<Integer> list = marker.getMarkedVertices();
		   if (list.size() == components.size()) {
			   return list;
		   }
		}
		return Collections.emptyList();
	}
	
	private void helixAlignment(List<List<Integer>> layerLines, int interactions) {
		if (layerLines.size() == 0 || layerLines.get(0).size() < 2) {
			return;
		}

		// get helix coordinates of subunit centroids
	//	Point3d[][] hl = subunitPairs(layerLines);
		
		// if there are less than 3 subunit centroids for superpostion,
		// use entire Calpha traces to superpose helix coordinates
	//	if (hl[0].length < 3) {
		Point3d[][]	hl = tracePairs(layerLines);
	//	}

		Matrix4d transformation = SuperPosition.superposeWithTranslation(hl[0], hl[1]);
		double rmsd = SuperPosition.rmsd(hl[0], hl[1]);
		System.out.println("HelixCheck.helixAlignment() rmsd: " + rmsd);

		if (rmsd < this.quatSymmetryParameters.getRmsdThreshold()) {
			HelixParameters hp = new HelixParameters(transformation, rmsd, interactions, layerLines);
			hp.setAlignmetMatrix(alignHelixAxis(hp));
//			System.out.println("Start count: " + getStartCount(hp));
			this.helixList.add(hp);
			System.out.println("Helix: " + this.helixList.get(this.helixList.size()-1));
		} 
	}
	
	private Point3d[][] tracePairs(List<List<Integer>> layerLines) {		
		List<Point3d[]> traces = this.subunits.getTraces();
		List<Point3d> xList = new ArrayList<Point3d>();
		List<Point3d> yList = new ArrayList<Point3d>();
		
		for (List<Integer> layerLine: layerLines) {					
//		for (int i = 0; i < layerLine.size()-1; i++) {
			for (int i = 0; i < 1; i++) {
			    System.out.println("HelixCheck.tracePairs: " + layerLine.get(i) + "-" + layerLine.get(i+1));
				for (Point3d p: traces.get(layerLine.get(i))) {
					xList.add(new Point3d(p));
				}
				for (Point3d p: traces.get(layerLine.get(i+1))) {
					yList.add(new Point3d(p));
				}
			}
			break;
		}
		Point3d[] x = xList.toArray(new Point3d[xList.size()]);
		Point3d[] y = yList.toArray(new Point3d[yList.size()]);
		Point3d[][] hl = {x, y};
		return hl;
	}
	
	private HelixParameters getByLowestTwistAngle() {
		double minAngle = Double.MAX_VALUE;
		HelixParameters p1 = null;
		
		for (HelixParameters p: this.helixList) {
			double angle = Math.abs(p.getAngle());
			if (angle < minAngle) {
				minAngle = angle;
				p1 = p;
			}
		}
		return p1;
	}
	
	private HelixParameters getByLargestInteraction() {
		int maxInteractions = 0;
		HelixParameters p1 = null;
		
		for (HelixParameters p: this.helixList) {
			int interactions = p.getInteractions();
			if (interactions > maxInteractions) {
				maxInteractions = interactions;
				p1 = p;
			}
		}
		return p1;
	}
	
	private HelixParameters getByLowestStartNumber() {
		int minStart = Integer.MAX_VALUE;
		HelixParameters p1 = null;
		
		for (HelixParameters p: this.helixList) {
			int start = p.getStart();
			if (start < minStart) {
				minStart = start;
				p1 = p;
			}
		}
		return p1;
	}

	private Matrix4d alignHelixAxis(HelixParameters helixParameters) {
		Vector3d[] axisVectors = new Vector3d[2];
		axisVectors[0] = helixParameters.getAxisVector();
		
		List<Point3d> centers = this.subunits.getOriginalCenters();
		Point3d centroid = this.subunits.getCentroid();
		
		int referenceSubunit = helixParameters.getLayerLines().get(0).get(0);
		Vector3d perpVector = new Vector3d();
		perpVector.sub(centroid, centers.get(referenceSubunit));
		perpVector.normalize();
		perpVector.cross(perpVector, axisVectors[0]);
		perpVector.normalize();
		perpVector.cross(perpVector, axisVectors[0]);
		perpVector.normalize();
		axisVectors[1] = perpVector;
		
		Vector3d[] referenceVectors = {Y_AXIS, Z_AXIS};    
		Matrix4d transformation = alignAxes(axisVectors, referenceVectors);
		
		// add translational component
//		Vector3d translation = new Vector3d(centroid);
//		translation.negate();
//		transformation.setTranslation(translation);
//		transformation.m33 = 1.0;
		return transformation;
	}
	
	private void sortLayerLines() {
		HelixParameters oneStart = getOneStart();
		if (oneStart == null) {
			return;
		}
		
		List<Integer> oneStartLine = oneStart.getLayerLines().get(0);
		
		for (HelixParameters p: this.helixList) {
			if (p.getStart() > 1) {
				List<List<Integer>> orderedList = new ArrayList<List<Integer>>();
				List<List<Integer>> layerLines = p.getLayerLines();
				System.out.println("original list: " + layerLines);
				for (int subunit: oneStartLine) {
					for (List<Integer> line: layerLines) {
						if (line.get(0) == subunit) {
							orderedList.add(line);
							break;
						}
					}
				}
				layerLines = orderedList;
				System.out.println("Ordered list: " + orderedList);
			}
		}
	}
	
	private HelixParameters getOneStart() {
		for (HelixParameters p: this.helixList) {
			if (p.getStart() == 1) {
				return p;
			}
		}
		return null;
	}
	
	private Point3d[] getLatticeVertices(Point3d[] points, Matrix4d transformation) {
		Matrix4d inverseTransformation = new Matrix4d();
		inverseTransformation.setIdentity();
		inverseTransformation.invert(transformation);
		
		Point3d[] latticePoints = new Point3d[points.length];
		for (int i = 0; i < points.length; i++) {
			latticePoints[i] = new Point3d(points[i]);
			transformation.transform(latticePoints[i]);
			double scale = (this.xzRadiusMax+8) / Math.sqrt(latticePoints[i].x*latticePoints[i].x + latticePoints[i].z*latticePoints[i].z);
     		latticePoints[i].x *= scale;
			latticePoints[i].z *= scale;
			inverseTransformation.transform(latticePoints[i]);
//			System.out.println("Lattice point: " + points[i] + " s: " + scale + " lp: " + latticePoints[i]);
		}
		return latticePoints;
	}
	
	private String getJmolLatticeLines(HelixParameters helixParameters, Matrix4d transformation, int index) {
		List<List<Integer>> layerLines = helixParameters.getLayerLines();
		List<Point3d> centers =  subunits.getOriginalCenters();
		Point3d pMin = new Point3d(minBoundary);
		Point3d pMax = new Point3d(maxBoundary);
		double width = pMin.distance(pMax)*0.005;

		StringBuilder s = new StringBuilder();
		for (List<Integer> line: layerLines) {
			s.append("draw lattice");
			s.append(index++);
			s.append(" line");
			Point3d[] points = new Point3d[line.size()];
			for (int i = 0; i < line.size(); i++) {
				points[i] = new Point3d(centers.get(line.get(i)));
			}
			points = getLatticeVertices(points, transformation);
			for (Point3d p: points) {
				s.append(getJmolPoint(p));
			}
			s.append("width ");
			s.append(fDot2(width));
			s.append(" color yellow");
			s.append(" off;");
		}
		
		return s.toString();
	}
	
	
	/**
	 * Returns a lower precision floating point number for Jmol
	 * @param f
	 * @return
	 */
	private static float jMolFloat(double f) {
	//	if (Math.abs(f) < 1.0E-7) {
	//		return 0.0f;
	//	}
		return (float)f;

	}
	
	private static String getJmolPoint(Tuple3d point) {
		StringBuilder s = new StringBuilder();
		s.append("{");
		s.append(fDot2(point.x));
		s.append(",");
		s.append(fDot2(point.y));
		s.append(",");
		s.append(fDot2(point.z));
		s.append("}");
		return s.toString();
	}
	
	private static String getJmolColor(Color4f color) {
		StringBuilder s = new StringBuilder();
		s.append("{");
		s.append(f1Dot2(color.x));
		s.append(",");
		s.append(f1Dot2(color.y));
		s.append(",");
		s.append(f1Dot2(color.z));
		s.append("}");
		return s.toString();
	}
	
	private void calcBoundaries(Matrix4d transformation) {	
		minBoundary.x = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.y = Double.MAX_VALUE;
		maxBoundary.x = Double.MIN_VALUE;
		minBoundary.z = Double.MAX_VALUE;
		maxBoundary.z = Double.MIN_VALUE;

		Point3d probe = new Point3d();

		for (Point3d[] list: subunits.getTraces()) {
			for (Point3d p: list) {
				probe.set(p);
				transformation.transform(probe);

				minBoundary.x = Math.min(minBoundary.x, probe.x);
				maxBoundary.x = Math.max(maxBoundary.x, probe.x);
				minBoundary.y = Math.min(minBoundary.y, probe.y);
				maxBoundary.y = Math.max(maxBoundary.y, probe.y);
				minBoundary.z = Math.min(minBoundary.z, probe.z);
				maxBoundary.z = Math.max(maxBoundary.z, probe.z);
				xzRadiusMax = Math.max(xzRadiusMax, Math.sqrt(probe.x*probe.x + probe.z*probe.z));
			}
		}
//		System.out.println("minBoundary: " + minBoundary);
//		System.out.println("maxBoundary: " + maxBoundary);
	}
	
	
	

	
	private boolean precheck () {
		if (! this.rotationGroup.getPointGroup().startsWith("C")) {
			return false;
		}
		int n = this.subunits.getSubunitCount();
		
		if (n < 3) {
			return false;
		}
		
		List<Integer> folds = this.subunits.getFolds();
		int maxFold = folds.get(folds.size()-1);
		System.out.println("maxFold : " + maxFold);

		return n % maxFold == 0;
	}
	
	
	/**
	 * Returns a transformation matrix that rotates refPoints to match
	 * coordPoints
	 * @param refPoints the points to be aligned
	 * @param referenceVectors
	 * @return
	 */
	private Matrix4d alignAxes(Vector3d[] axisVectors, Vector3d[] referenceVectors) {
		Matrix4d m1 = new Matrix4d();
		AxisAngle4d a = new AxisAngle4d();
		Vector3d axis = new Vector3d();
		
		// calculate rotation matrix to rotate refPoints[0] into coordPoints[0]
		Vector3d v1 = new Vector3d(axisVectors[0]);
		Vector3d v2 = new Vector3d(referenceVectors[0]);
		double dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m1.set(a);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m1.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around x-axis
			m1.set(flipX());
		}
		
		// apply transformation matrix to all refPoints
		m1.transform(axisVectors[0]);
		m1.transform(axisVectors[1]);
		
		// calculate rotation matrix to rotate refPoints[1] into coordPoints[1]
		v1 = new Vector3d(axisVectors[1]);
		v2 = new Vector3d(referenceVectors[1]);
		Matrix4d m2 = new Matrix4d();
		dot = v1.dot(v2);
		if (Math.abs(dot) < 0.999) {
			axis.cross(v1,v2);
			axis.normalize();
			a.set(axis, v1.angle(v2));
			m2.set(a);
		} else if (dot > 0) {
			// parallel axis, nothing to do -> identity matrix
			m2.setIdentity();
		} else if (dot < 0) {
			// anti-parallel axis, flip around z-axis
			m2.set(flipZ());
		}
		
		// apply transformation matrix to all refPoints
		m2.transform(axisVectors[0]);
		m2.transform(axisVectors[1]);
		
		// combine the two rotation matrices
		m2.mul(m1);

		// the RMSD should be close to zero
		Point3d[] axes = new Point3d[2];
		axes[0] = new Point3d(axisVectors[0]);
		axes[1] = new Point3d(axisVectors[1]);
		Point3d[] ref = new Point3d[2];
		ref[0] = new Point3d(referenceVectors[0]);
		ref[1] = new Point3d(referenceVectors[1]);
		if (SuperPosition.rmsd(axes, ref) > 0.01) {
			System.out.println("Warning: AxisTransformation: axes alignment is off. RMSD: " + SuperPosition.rmsd(axes, ref));
		}
		
		return m2;
	}
	
	private static String f1Dot2(float number) {
		return String.format("%1.2f", number);
	}
	
	private static String fDot2(double number) {
		return String.format("%.2f", number);
	}
	
	private static Matrix4d flipX() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = 1;
		rot.m11 = -1;
		rot.m22 = -1;
		rot.m33 = 1;
		return rot;
	}
	
	private static Matrix4d flipZ() {
		Matrix4d rot = new Matrix4d();
		rot.m00 = -1;
		rot.m11 = -1;
		rot.m22 = 1;
		rot.m33 = 1;
		return rot;
	}
}
