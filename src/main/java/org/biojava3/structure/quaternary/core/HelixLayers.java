/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Peter
 */
public class HelixLayers {
	private List<Helix> helices = new ArrayList<Helix>();

	public int size() {
		return helices.size();
	}

	public void addHelix(Helix helix) {
		helices.add(helix);
	}
	
	public Helix getHelix(int index) {
		return helices.get(index);
	}
	
	/*
	 * Returns Helix with lowest twist angle
	 */
	public Helix getByLowestAngle() {
		double angle = Double.MAX_VALUE;
		Helix lowest = null;
		for (Helix helix: helices) {
			if (helix.getAngle() < angle) {
				angle = helix.getAngle();
				lowest = helix;
			}
		}
		return lowest;
	}
	
	/*
	 * Returns Helix with largest number of intermolecular contacts
	 * between repeat units
	 */
	public Helix getByLargestContacts() {
		double contacts = 0;
		Helix largest = null;
		for (Helix helix: helices) {
			if (helix.getContacts() > contacts) {
				contacts = helix.getContacts();
				largest = helix;
			}
		}
		return largest;
	}
	
	/* 
	 * Returns Helix that has the largest number of contacts, besides
	 * the Helix with the lowest twist angle
	 */
	public Helix getByLargestContactsNotLowestAngle() {
		double contacts = 0;
		Helix lowest = getByLowestAngle();
		// TODO why are there helices with almost identical helix parameters??
		double angle = lowest.getAngle() + 0.05;
		Helix largest = null;
		for (Helix helix: helices) {
			if (helix == lowest) {
				continue;
			}
			if (helix.getContacts() > contacts && helix.getAngle() > angle) {
				contacts = helix.getContacts();
				largest = helix;
			}
		}
		if (largest == null) {
			return lowest;
		}
		return largest;
	}
	
	public double getAverageSubunitRmsd() {
		if (size() == 0) {
			return 0;
		}
		double rmsd = 0;
		for (Helix helix: helices) {
			rmsd+= helix.getSubunitRmsd();
		}
		return rmsd/this.size();
	}

	public double getAverageTraceRmsd() {
		if (size() == 0) {
			return 0;
		}
		double rmsd = 0;
		for (Helix helix: helices) {
			rmsd+= helix.getTraceRmsd();
		}
		return rmsd/this.size();
	}
	
	public double getAverageTraceTmScoreMin() {
		if (size() == 0) {
			return 0;
		}
		double sum = 0;
		for (Helix helix: helices) {
			sum+= helix.getTraceTmScoreMin();
		}
		return sum/this.size();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Helices: " + size() + "\n");
		for (Helix s: helices) {
			sb.append(s.toString() + "\n");
		}
		return sb.toString();
	}

}
