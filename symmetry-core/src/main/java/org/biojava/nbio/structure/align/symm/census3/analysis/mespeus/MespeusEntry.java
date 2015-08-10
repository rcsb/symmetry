package org.biojava.nbio.structure.align.symm.census3.analysis.mespeus;

/**
 * A metalloprotein from <a href="http://mespeus.bch.ed.ac.uk/MESPEUS">MESPEUS</a>.
 * @author dmyersturnbull
 */
public class MespeusEntry {

	private final int id;
	private final float distance;
	private final int coordinationNumber;
	private final CoordinationGeometry shape;
	private final String metalName;
	private final String donorResidueName;
	private final String pdbId;
	private final double resolution;
	private final float rmsd;
	private final int difference;
	public static MespeusEntry parse(String line) {
		String[] parts = line.split("\t");
		return new MespeusEntry(Integer.parseInt(parts[0]), Float.parseFloat(parts[1]), Integer.parseInt(parts[2]), CoordinationGeometry.parse(parts[3]), parts[4], parts[5], parts[6], Float.parseFloat(parts[7]), Float.parseFloat(parts[8]), Integer.parseInt(parts[9]));
	}
	public MespeusEntry(int id, float distance, int coordinationNumber, CoordinationGeometry shape, String metalName,
			String donorResidueName, String pdbId, float resolution, float rmsd, int difference) {
		super();
		this.id = id;
		this.distance = distance;
		this.coordinationNumber = coordinationNumber;
		this.shape = shape;
		this.metalName = metalName;
		this.donorResidueName = donorResidueName;
		this.pdbId = pdbId.toLowerCase();
		this.resolution = resolution;
		this.rmsd = rmsd;
		this.difference = difference;
	}
	public int getId() {
		return id;
	}
	public float getDistance() {
		return distance;
	}
	public int getCoordinationNumber() {
		return coordinationNumber;
	}
	public CoordinationGeometry getShape() {
		return shape;
	}  
	public String getMetalName() {
		return metalName;
	}
	public String getDonorResidueName() {
		return donorResidueName;
	}
	public String getPdbId() {
		return pdbId;
	}
	public double getResolution() {
		return resolution;
	}
	public float getRmsd() {
		return rmsd;
	}
	public int getDifference() {
		return difference;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + coordinationNumber;
		result = prime * result + difference;
		result = prime * result + Float.floatToIntBits(distance);
		result = prime * result + ((donorResidueName == null) ? 0 : donorResidueName.hashCode());
		result = prime * result + id;
		result = prime * result + ((metalName == null) ? 0 : metalName.hashCode());
		result = prime * result + ((pdbId == null) ? 0 : pdbId.hashCode());
		long temp;
		temp = Double.doubleToLongBits(resolution);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + Float.floatToIntBits(rmsd);
		result = prime * result + ((shape == null) ? 0 : shape.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		MespeusEntry other = (MespeusEntry) obj;
		if (coordinationNumber != other.coordinationNumber) return false;
		if (difference != other.difference) return false;
		if (Float.floatToIntBits(distance) != Float.floatToIntBits(other.distance)) return false;
		if (donorResidueName == null) {
			if (other.donorResidueName != null) return false;
		} else if (!donorResidueName.equals(other.donorResidueName)) return false;
		if (id != other.id) return false;
		if (metalName == null) {
			if (other.metalName != null) return false;
		} else if (!metalName.equals(other.metalName)) return false;
		if (pdbId == null) {
			if (other.pdbId != null) return false;
		} else if (!pdbId.equals(other.pdbId)) return false;
		if (Double.doubleToLongBits(resolution) != Double.doubleToLongBits(other.resolution)) return false;
		if (Float.floatToIntBits(rmsd) != Float.floatToIntBits(other.rmsd)) return false;
		if (shape == null) {
			if (other.shape != null) return false;
		} else if (!shape.equals(other.shape)) return false;
		return true;
	}
	@Override
	public String toString() {
		return "MespeusEntry [id=" + id + ", distance=" + distance + ", coordinationNumber=" + coordinationNumber
				+ ", shape=" + shape + ", metalName=" + metalName + ", donorResidueName=" + donorResidueName
				+ ", pdbId=" + pdbId + ", resolution=" + resolution + ", rmsd=" + rmsd + ", difference=" + difference
				+ "]";
	}
	
}