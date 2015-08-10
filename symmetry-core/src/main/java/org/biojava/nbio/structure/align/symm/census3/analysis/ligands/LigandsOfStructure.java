package org.biojava.nbio.structure.align.symm.census3.analysis.ligands;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import javax.xml.bind.annotation.XmlTransient;

/**
 * A list of ligands associated with a single structure.
 * 
 * @author dmyersturnbull
 */
public class LigandsOfStructure implements Serializable {

	private static final long serialVersionUID = 8844385971743419028L;
	private List<CensusLigand> ligands = new ArrayList<CensusLigand>();
	private String structureName;

	public LigandsOfStructure() {
		super();
	}

	public LigandsOfStructure(String structureName) {
		super();
		this.structureName = structureName;
	}

	public boolean add(CensusLigand e) {
		return ligands.add(e);
	}

	public boolean addAll(Collection<? extends CensusLigand> c) {
		return ligands.addAll(c);
	}

	public boolean contains(Object o) {
		return ligands.contains(o);
	}

	public CensusLigand get(int index) {
		return ligands.get(index);
	}

	public List<CensusLigand> getLigands() {
		return ligands;
	}

	@XmlTransient
	public String getStructureName() {
		return structureName;
	}

	public boolean isEmpty() {
		return ligands.isEmpty();
	}

	public Iterator<CensusLigand> iterator() {
		return ligands.iterator();
	}

	public void setLigands(List<CensusLigand> ligands) {
		this.ligands = ligands;
	}

	public void setStructureName(String structureName) {
		this.structureName = structureName;
	}

	public int size() {
		return ligands.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		sb.append(structureName);
		if (!ligands.isEmpty())
			sb.append(": ");
		for (int i = 0; i < ligands.size(); i++) {
			sb.append(ligands.get(i));
			if (i < ligands.size() - 1)
				sb.append(", ");
		}
		sb.append("]");
		return sb.toString();
	}

}
