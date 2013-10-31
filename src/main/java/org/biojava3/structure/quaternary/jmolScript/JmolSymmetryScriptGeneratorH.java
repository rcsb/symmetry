/**
 * 
 */
package org.biojava3.structure.quaternary.jmolScript;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Color4f;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

import org.biojava3.structure.quaternary.core.HelixAxisAligner;
import org.biojava3.structure.quaternary.core.Subunits;
import org.biojava3.structure.quaternary.utils.ColorBrewer;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorH extends JmolSymmetryScriptGenerator { 
	private static double AXIS_SCALE_FACTOR = 1.2;
	private static double SIDE_CHAIN_EXTENSION = 6.0;
	private HelixAxisAligner helixAxisAligner = null;
	private String name = "";
	private String defaultColoring = "";

	public JmolSymmetryScriptGeneratorH(HelixAxisAligner helixAxisAligner, String name) {
        this.helixAxisAligner = helixAxisAligner;
        this.name = name;
	}
	
	public int getZoom() {
		int zoom = (int) Math.round(100.0/AXIS_SCALE_FACTOR);
		return zoom;
	}
	
	/**
	 * Returns a Jmol script to set the default orientation for a structure
	 * @return Jmol script
	 */
	public String getDefaultOrientation() {	
		StringBuilder s = new StringBuilder();
		s.append(setCentroid());

		Quat4d q = new Quat4d();
		q.set(helixAxisAligner.getRotationMatrix());
		
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
	
	public int getOrientationCount() {
		return 4;
	}
	
	@Override
	public String getOrientation(int index) {
		StringBuilder s = new StringBuilder();
		s.append(setCentroid());

		// calculate  orientation
		Matrix4d matrix = new Matrix4d();
		matrix.setIdentity();
		matrix.setRotation(new AxisAngle4d(1,0,0,Math.PI/2*(index+1)));
		Quat4d r = new Quat4d();
		r.set(new AxisAngle4d(1,0,0,index*Math.PI/2));
		Quat4d q = new Quat4d();
		q.set(helixAxisAligner.getRotationMatrix());	
		r.mul(q);

		// set orientation
		s.append("moveto 4 quaternion{");
		s.append(jMolFloat(r.x));
		s.append(",");
		s.append(jMolFloat(r.y));
		s.append(",");
		s.append(jMolFloat(r.z));
		s.append(",");
		s.append(jMolFloat(r.w));
		s.append("}");
		s.append(";");
		return s.toString();
	}
	
	/**
	 * Returns a Jmol script that sets a specific orientation and zoom
	 * to draw either axes or polyhedron
	 * @param index orientation index
	 * @return Jmol script
	 */
	public String getOrientationWithZoom(int index) {
		StringBuilder s = new StringBuilder();
		s.append(getOrientation(index));
		s.insert(s.length()-1, getZoom());
		return s.toString();
		
	}
	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	public String getOrientationName(int index) {	
		switch (index) {
		    case 0: return "Side";
		    case 1: return "Front";
		    case 2: return "Other side";
		    case 3: return "Back";
		    default: return "";
		}
	}

	/* (non-Javadoc)
	 * @see org.biojava3.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getTransformation()
	 */
	public Matrix4d getTransformation() {	
	    return helixAxisAligner.getTransformation();
	}
	
	
	public void setDefaultColoring(String colorScript) {
		this.defaultColoring = colorScript;
	}
	
	/**
	 * Returns a Jmol script that draws an invisible polyhedron around a structure.
	 * Use showPolyhedron() and hidePolyhedron() to toggle visibility.
	 * @return Jmol script
	 */
	public String drawPolyhedron() {
		StringBuilder s = new StringBuilder();
		
		// for now, return empty script until this has been implemented properly
//		if (s.length() == 0)
//			return s.toString();

//		Point3d[] vertices = getRepeatUnitCenters();
		List<Point3d> vertices = helixAxisAligner.getSubunits().getOriginalCenters();
		vertices = extendUnitCenters(vertices);
		
		int index = 0;
		Color4f c = new Color4f(Color.MAGENTA);
//		double width = getMaxExtension()*0.015;
		double width = getMaxExtension()*0.007;

//		List<List<Integer>> layerLines = helixAxisAligner.getHelixLayers().getByLargestContactsNotLowestAngle().getLayerLines();
//		layerLines.addAll(helixAxisAligner.getHelixLayers().getByLowestAngle().getLayerLines());
		List<List<Integer>> layerLines = helixAxisAligner.getHelixLayers().getByLargestContacts().getLayerLines();

		
		for (List<Integer> line : layerLines) {
			s.append("draw polyhedron");
			s.append(name);
			s.append(index++);
			s.append(" line");
			for (int i: line) {
				s.append(getJmolPoint(vertices.get(i)));
			}
			s.append("width ");
		    s.append(fDot2(width));
			s.append(" color");
			s.append(getJmolColor(c));
			s.append(" off;");
		}
		
		List<Point3d> interiorVertices = interiorCenters(vertices);
		
		for (int i = 0; i < vertices.size(); i++) {
			s.append("draw polyhedron");
			s.append(name);
			s.append(index++);
			s.append(" line");
			s.append(getJmolPoint(vertices.get(i)));
			s.append(getJmolPoint(interiorVertices.get(i)));
			s.append("width ");
		    s.append(fDot2(width));
			s.append(" color");
			s.append(getJmolColor(c));
			s.append(" off;");
		}
		
//		Matrix4d transformation = helixAxisTransformation.getHelixLayers().getByLowestAngle().getTransformation();
//		
//		Point3d[] vertices2 = SuperPosition.clonePoint3dArray(vertices);
//		for (Point3d p: vertices2) {
//			transformation.transform(p);
//		}
//		
//		// extend grid by one turn in +z direction
//		for (int[] lineLoop: getLayerLines()) {
//			s.append("draw line");
//			s.append(name);
//			s.append(index++);
//			s.append(" line");
//			int count = 0;
//			for (int i: lineLoop) {
//				for (Point3d v: vertices) {
//					if (v.distance(vertices2[i]) > 0.1) {
//						count++;
//						break;
//					}
//				}
//				if (count > 0) {
//					s.append(getJmolPoint(vertices2[i]));
//				} else {
//					s.append(" ");
//				}
//			}
//			s.append("width ");
//		    s.append(fDot2(width));
//			s.append(" color");
//			s.append(getJmolColor(c));
//			s.append(" off;");
//		}
//		
//		
//		// extend grid by one turn into -z direction
//		Point3d[] vertices3 = SuperPosition.clonePoint3dArray(vertices);
//		Matrix4d m = new Matrix4d();
//		m.set(transformation);
//		m.invert();
//		for (Point3d p: vertices3) {
//			m.transform(p);
//		}
//
//		for (int[] lineLoop: getLayerLines()) {
//			s.append("draw line");
//			s.append(name);
//			s.append(index++);
//			s.append(" line");
//			int count = 0;
//			for (int i: lineLoop) {
//				for (Point3d v: vertices) {
//					if (v.distance(vertices3[i]) > 0.1) {
//						count++;
//					}
//				}
//				if (count > 0) {
//					s.append(getJmolPoint(vertices3[i]));
//				} else {
//					s.append(" ");
//				}
//			}
//			s.append("width ");
//		    s.append(fDot2(width));
//			s.append(" color");
//			s.append(getJmolColor(c));
//			s.append(" off;");
//		}
		
		return s.toString();
	}
	
	public String hidePolyhedron() {
		return "draw polyhedron* off;";
	}
	
	public String showPolyhedron() {
		return "draw polyhedron* on;";
	}
	
	/**
	 * Returns a Jmol script that draws symmetry or inertia axes for a structure.
	 * Use showAxes() and hideAxes() to toggle visibility.
	 * @return Jmol script
	 */
	public String drawAxes() {
		StringBuilder s = new StringBuilder();
//		Point3d centroid = helixAxisAligner.getCentroid();
//		System.out.println("Centroid: " + helixAxisAligner.getCentroid());
//		Point3d centroid = helixAxisAligner.getGeometricCenter();
		Point3d centroid = helixAxisAligner.calcCenterOfRotation();
//		System.out.println("Geometric center: " + centroid);
		AxisAngle4d axisAngle = helixAxisAligner.getHelixLayers().getByLowestAngle().getAxisAngle();
		Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);

		s.append("draw axesHelical");
		s.append(name);
		s.append(0);
		s.append(" ");
		s.append("line");
		Point3d v1 = new Point3d(axis);
		v1.scale(AXIS_SCALE_FACTOR*(helixAxisAligner.getDimension().y + SIDE_CHAIN_EXTENSION));
		Point3d v2 = new Point3d(v1);
		v2.negate();
		v1.add(centroid);
		v2.add(centroid);
		s.append(getJmolPoint(v1));
		s.append(getJmolPoint(v2));
		s.append("width 1.0 ");
		s.append(" color red");
		s.append(" off;");
		return s.toString();
	}
	
	/**
	 * Returns a Jmol script to hide axes
	 * @return Jmol script
	 */
	public String hideAxes() {
		return "draw axes* off;";
	}
	
	/**
	 * Returns a Jmol script to show axes
	 * @return Jmol script
	 */
	public String showAxes() {
		return "draw axes* on;";
	}
	
	/**
	 * Returns a Jmol script that displays a symmetry polyhedron and symmetry axes
	 * and then loop through different orientations
	 * @return Jmol script
	 */
	public String playOrientations() {
		StringBuilder s = new StringBuilder();
		
		// draw footer	
		s.append(drawFooter("Symmetry Helical", "white"));
		
		// draw polygon
		s.append(drawPolyhedron()); // draw invisibly
		s.append(showPolyhedron());
			
		// draw axes
		s.append(drawAxes());
		s.append(showAxes());
		
		// loop over all orientations with 4 sec. delay
		for (int i = 0; i < getOrientationCount(); i++) {
			s.append(deleteHeader());
			s.append(getOrientationWithZoom(i));
			s.append(drawHeader(getOrientationName(i), "white"));
			s.append("delay 4;");
		}
		
		// go back to first orientation
		s.append(deleteHeader());
		s.append(getOrientationWithZoom(0));
		s.append(drawHeader(getOrientationName(0), "white"));
		
		return s.toString();
	}
	
	
	/**
	 * Returns a Jmol script that colors the subunits of a structure by different colors
	 * @return
	 */
	public String colorBySubunit() {
	    Subunits subunits = helixAxisAligner.getSubunits();
	    List<Integer> modelNumbers = subunits.getModelNumbers();
	    List<String> chainIds = subunits.getChainIds();
	    List<List<Integer>> orbits = helixAxisAligner.getOrbits();

		Color4f[] colors = ColorBrewer.Spectral.getColor4fPalette(orbits.size());
		int half = colors.length/2;
		for (int i = 0; i < half; i++) {
			if (i % 2 != 0) {
			   Color4f temp = colors[i];
			   colors[i] = colors[half+i];
			   colors[half+i] = temp;
			}
		}
	    Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
	    
		for (int i = 0; i < orbits.size(); i++) {
				for (Integer j: orbits.get(i)) {
					Color4f c = colors[i];
					List<String> ids = colorMap.get(c);
					if (ids == null) {
						ids = new ArrayList<String>();
						colorMap.put(c,  ids);
					}
					String id = chainIds.get(j) + "/" + (modelNumbers.get(j)+1);
					ids.add(id);
				}
		}
		String coloring = defaultColoring + getJmolColorScript(colorMap);
		return coloring;
	}
	
	/**
	 * Returns a Jmol script that colors subunits by their sequence cluster ids.
	 * @return Jmol script
	 */
	public String colorBySequenceCluster() {
	    Subunits subunits = helixAxisAligner.getSubunits();
	    int n = subunits.getSubunitCount();
	    List<Integer> modelNumbers = subunits.getModelNumbers();
	    List<String> chainIds = subunits.getChainIds();
	    List<Integer> seqClusterIds = subunits.getSequenceClusterIds();
	    int clusters = Collections.max(seqClusterIds) + 1;
	    Color4f[] colors = ColorBrewer.BrBG.getColor4fPalette(clusters);
		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
		
		for (int i = 0; i < n; i++) {
			Color4f c = colors[seqClusterIds.get(i)];
			List<String> ids = colorMap.get(c);
			if (ids == null) {
				ids = new ArrayList<String>();
				colorMap.put(c,  ids);
			}
			String id = chainIds.get(i) + "/" + (modelNumbers.get(i)+1);
			ids.add(id);

		}
		String coloring = defaultColoring + getJmolColorScript(colorMap);
		return coloring;
	}
	
	
	/**
	 * Returns a Jmol script that colors subunits to highlight the symmetry within a structure
	 * Different subunits should have a consistent color scheme or different shade of the same colors
	 * @return Jmol script
	 */
	public String colorBySymmetry() {
		List<List<Integer>> units = helixAxisAligner.getHelixLayers().getByLargestContacts().getLayerLines();
		units = orientLayerLines(units);
		Subunits subunits = helixAxisAligner.getSubunits();
		List<Integer> modelNumbers = subunits.getModelNumbers();
		List<String> chainIds = subunits.getChainIds();
		List<Integer> clusterIds = subunits.getSequenceClusterIds();
		int clusterCount = Collections.max(clusterIds) + 1;

		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>(); 

		int maxLen = 0;
		for (List<Integer> unit: units) {
			maxLen = Math.max(maxLen,  unit.size());
		}
	
//		Color4f[] colors = getSymmetryColors(permutation.size()); 
		Color4f[] colors = getSymmetryColors(subunits.getSubunitCount()); 
		int count = 0;
		
		for (int i = 0; i < maxLen; i++) {
			for (int j = 0; j < units.size(); j++) {
				int m = units.get(j).size();
				if (i < m) {
					int subunit = units.get(j).get(i);
					int cluster = clusterIds.get(subunit);
					float scale = 0.3f + 0.7f * (float) (cluster+1)/clusterCount;
					Color4f c = new Color4f(colors[count]);
					count++;
					c.scale(scale);
					List<String> ids = colorMap.get(c);
					if (ids == null) {
						ids = new ArrayList<String>();
						colorMap.put(c,  ids);
					}
					String id = chainIds.get(subunit) + "/" + (modelNumbers.get(subunit)+1);
					ids.add(id);
				}
			}
		}
		
		String coloring = defaultColoring + getJmolColorScript(colorMap);
		return coloring;
	}
	
	/**
	 * Orients layer lines from lowest y-axis value to largest y-axis value
	 */
	private List<List<Integer>> orientLayerLines(List<List<Integer>> layerLines) {
		Matrix4d transformation = helixAxisAligner.getTransformation();
		List<Point3d> centers = helixAxisAligner.getSubunits().getOriginalCenters();
		
		for (int i = 0; i < layerLines.size(); i++) {
			List<Integer> layerLine = layerLines.get(i);
			
			// get center of first subunit in layerline and transform to standard orientation (helix axis aligned with y-axis)
			int first = layerLine.get(0);
			Point3d firstSubunit = new Point3d(centers.get(first));
			transformation.transform(firstSubunit);
			
			// get center of last subunit in layerline and transform to standard orientation (helix axis aligned with y-axis)
			int last = layerLine.get(layerLine.size()-1);
			Point3d lastSubunit = new Point3d(centers.get(last));
			transformation.transform(lastSubunit);
			
			// a layerline should start at the lowest y-value, so all layerlines have a consistent direction from -y value to +y value
			if (firstSubunit.y > lastSubunit.y) {
//				System.out.println("reorienting layer line: " + layerLine);
				Collections.reverse(layerLine);
			}
		}
		return layerLines;
	}
		
	// --- protected methods ---
	/**
	 * Returns the maximum extension (length) of structure
	 * @return
	 */
	protected double getMaxExtension() {
		Vector3d dimension = helixAxisAligner.getDimension();
		double maxExtension = Math.max(dimension.x, dimension.y);
		maxExtension = Math.max(maxExtension, dimension.z);
		return maxExtension;
	}
	
	private List<Point3d> extendUnitCenters(List<Point3d> points) {
		Matrix4d transformation = helixAxisAligner.getTransformation();
		Matrix4d reversetransformation = helixAxisAligner.getReverseTransformation();
		double radius = helixAxisAligner.getRadius();
		radius += SIDE_CHAIN_EXTENSION; // add space for side chains
//		System.out.println("XZRadius: " + radius);
		

		List<Point3d> ePoints = new ArrayList<Point3d>(points.size());
		for (Point3d p: points) {
			Point3d q = new Point3d(p);
            transformation.transform(q);
            double d = Math.sqrt(q.x*q.x + q.z*q.z);
            double scale = radius/d;
            q.x *= scale;
            q.z *= scale;
            reversetransformation.transform(q);
            ePoints.add(q);
		}
		return ePoints;
	}
	
	private List<Point3d> interiorCenters(List<Point3d> points) {
		Matrix4d transformation = helixAxisAligner.getTransformation();
		Matrix4d reversetransformation = helixAxisAligner.getReverseTransformation();
		
		List<Point3d> ePoints = new ArrayList<Point3d>(points.size());
		for (Point3d p: points) {
			Point3d q = new Point3d(p);
            transformation.transform(q);
            q.x = 0;
            q.z = 0;
            reversetransformation.transform(q);
            ePoints.add(q);
		}
		return ePoints;
	}
	
	private String drawHeader(String text, String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo top center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo ");
		s.append(text);
		s.append(";");
		return s.toString();
	}
	
	private String deleteHeader() {
		return "set echo top center;echo ;";
	}
	
	private String drawFooter(String text, String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo bottom center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo "+ text);
		s.append(";");
		return s.toString();
	}
	
	private String getJmolColorScript(Map<Color4f, List<String>> map) {
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
	 * Returns a unique color palette based on point group
	 * @param nColors
	 * @return
	 */
	private Color4f[] getSymmetryColors(int nColors) {
		int offset = 0;
		int dMax = nColors + offset;
		Color4f[] colors = ColorBrewer.Spectral.getColor4fPalette(dMax);
		System.arraycopy(colors, offset, colors, 0, dMax-offset);
		return colors;	
	}
	
 	private String setCentroid() {
		// calculate center of rotation
	//	Point3d centroid = axisTransformation.getGeometricCenter();
	//	Point3d centroid = helixAxisAligner.getCentroid();
		Point3d centroid = helixAxisAligner.calcCenterOfRotation();
			
		// set centroid
		StringBuilder s = new StringBuilder();
		s.append("center");
		s.append(getJmolPoint(centroid));
		s.append(";");
		return s.toString();
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
	
	private static String f1Dot2(float number) {
		return String.format("%1.2f", number);
	}
	
	private static String fDot2(double number) {
		return String.format("%.2f", number);
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

	
	
	
}
