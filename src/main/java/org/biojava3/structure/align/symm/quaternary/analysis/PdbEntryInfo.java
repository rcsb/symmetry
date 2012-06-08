package org.biojava3.structure.align.symm.quaternary.analysis;

import java.util.Date;
import java.util.List;

public class PdbEntryInfo {
	public enum ExperimentalMethod {ELECTRON_CRYSTALLOGRAPHY, ELECTRON_MICROSCOPY, FIBER_DIFFRACTION, NEUTRON_DIFFRACTION, SOLID_STATE_NMR, SOLUTION_NMR, SOLUTION_SCATTERING, THEORETICAL_MODEL, X_RAY,};
	private String pdbId;
	private int bioAssemblyCount;
	private Date releaseDate;
	private float resolution;
	private List<ExperimentalMethod> experimentalMethods;
	private List<String> entityTypes;
	
	/**
	 * @return the pDBID
	 */
	public String getPdbId() {
		return pdbId;
	}
	/**
	 * @param pDBID the pDBID to set
	 */
	public void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}
	/**
	 * @return the bioAssemblyCount
	 */
	public int getBioAssemblyCount() {
		return bioAssemblyCount;
	}
	/**
	 * @param bioAssemblyCount the bioAssemblyCount to set
	 */
	public void setBioAssemblyCount(int bioAssemblyCount) {
		this.bioAssemblyCount = bioAssemblyCount;
	}
	/**
	 * @return the releaseDate
	 */
	public Date getReleaseDate() {
		return releaseDate;
	}
	/**
	 * @param releaseDate the releaseDate to set
	 */
	public void setReleaseDate(Date releaseDate) {
		this.releaseDate = releaseDate;
	}
	/**
	 * @return the resolution
	 */
	public float getResolution() {
		return resolution;
	}
	/**
	 * @param resolution the resolution to set
	 */
	public void setResolution(float resolution) {
		this.resolution = resolution;
	}
	/**
	 * @return the experimentalMethods
	 */
	public List<ExperimentalMethod> getExperimentalMethods() {
		return experimentalMethods;
	}
	/**
	 * @param experimentalMethods the experimentalMethods to set
	 */
	public void setExperimentalMethods(List<ExperimentalMethod> experimentalMethods) {
		this.experimentalMethods = experimentalMethods;
	}
	/**
	 * @return the entityTypes
	 */
	public List<String> getEntityTypes() {
		return entityTypes;
	}
	/**
	 * @param entityTypes the entityTypes to set
	 */
	public void setEntityTypes(List<String> entityTypes) {
		this.entityTypes = entityTypes;
	}
	
	public boolean isProtein() {
		for (String s: entityTypes) {
			if (! s.equals("protein")) 
				return false;
		}
		return true;
	}
	
	public boolean containsProtein() {
		for (String s: entityTypes) {
			if (s.equals("protein")) 
				return true;
		}
		return false;
	}
	
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append(pdbId);
		s.append(" bioAssemblies: ");
		s.append(bioAssemblyCount);
		s.append(" release date: ");
		s.append(releaseDate);
		s.append(" resolution: ");
		s.append(resolution);
		s.append(" experimental method: ");
		s.append(experimentalMethods);
		s.append(" entity types: ");
		s.append(entityTypes);
		return s.toString();
	}
}
