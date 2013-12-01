package org.biojava3.structure.align.symm.census2.analysis;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The shape of a coordination complex. Contains deviations from ideal angles of ideal coordination complexes of the same coordination number.
 * @author dmyersturnbull
 */
public class CoordinationGeometry {

	private Map<CoordinationGeometryType, Float> deltas;
	
	private static final String REGEX = "(?:δ[ ]*([a-z]+)\\(([0-9.]+)°\\)).*?(?:δ[ ]*([a-z]+)\\(([0-9.]+)°\\))?";

	public static CoordinationGeometry parse(String field) {
		CoordinationGeometry shape = new CoordinationGeometry();
		shape.deltas = new LinkedHashMap<CoordinationGeometryType, Float>();
		Pattern pattern = Pattern.compile(REGEX);
		Matcher matcher = pattern.matcher(field);
		matcher.find();
		CoordinationGeometryType geometryType = null;
		if (matcher.matches()) {
			for (int i = 1; i <= matcher.groupCount(); i++) {
				if (i % 2 == 1) {
					String geometry = matcher.group(i);
					if (geometry == null) break;
					if (geometry.equals("oct")) {
						geometryType = CoordinationGeometryType.OCTAHEDRAL;
					} else if (geometry.equals("tet")) {
						geometryType = CoordinationGeometryType.TETRAHEDRAL;
					} else if (geometry.equals("tetp")) {
						geometryType = CoordinationGeometryType.SQUARE_PYRAMIDAL;
					} else if (geometry.equals("sqp")) {
						geometryType = CoordinationGeometryType.SQUARE_PLANAR;
					} else if (geometry.equals("tbp")) {
						geometryType = CoordinationGeometryType.TRIGONAL_BIPYRAMIDAL;
					} else {
						throw new IllegalArgumentException("Couldn't understand geometry " + geometry);
					}
				} else {
					float delta = Float.parseFloat(matcher.group(i));
					shape.deltas.put(geometryType, delta);
					geometryType = null; // just for safety
				}
			}
		}
		return shape;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		// deltas.hashCode() checks the entrySet, so we're fine
		result = prime * result + ((deltas == null) ? 0 : deltas.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		CoordinationGeometry other = (CoordinationGeometry) obj;
		if (deltas == null) {
			if (other.deltas != null) return false;
		}
		if (deltas.size() != other.deltas.size()) {
			return false;
		}
		for (Map.Entry<CoordinationGeometryType, Float> entry : deltas.entrySet()) {
			if (!other.deltas.containsKey(entry.getKey()) || other.deltas.get(entry.getKey()) != entry.getValue()) return false;
		}
		for (Map.Entry<CoordinationGeometryType, Float> entry : other.deltas.entrySet()) {
			if (!deltas.containsKey(entry.getKey()) || deltas.get(entry.getKey()) != entry.getValue()) return false;
		}
		return true;
	}

	public Float getDelta(CoordinationGeometryType key) {
		return deltas.get(key);
	}

	public Map<CoordinationGeometryType, Float> getDeltas() {
		return deltas;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("[");
		for (Map.Entry<CoordinationGeometryType, Float> entry : deltas.entrySet()) {
			sb.append("{δ(" + entry.getKey().name().toLowerCase() + ")=" + entry.getValue() + "°}");
		}
		sb.append("]");
		return sb.toString();
	}
	
}