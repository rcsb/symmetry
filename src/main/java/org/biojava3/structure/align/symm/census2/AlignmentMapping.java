package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.annotation.adapters.XmlAdapter;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;

/**
 * An alignment function.
 * In future versions, may be able to handle different mappings between different symmetry subunits.
 * @author dmyersturnbull
 */
public class AlignmentMapping implements Serializable {

	private static final long serialVersionUID = -7879437940997148046L;

	private Map<Integer, Integer> simpleFunction;

	public AlignmentMapping() {
	}

	public AlignmentMapping(AFPChain afpChain) {
		try {
			this.simpleFunction = AlignmentTools.alignmentAsMap(afpChain);
		} catch (StructureException e) {
			throw new RuntimeException(e);
		}
	}
	
	@XmlJavaTypeAdapter(AlignmentFunctionAdapter.class)
	public Map<Integer, Integer> getSimpleFunction() {
		return simpleFunction;
	}

	public void setSimpleFunction(Map<Integer, Integer> simpleFunction) {
		this.simpleFunction = simpleFunction;
	}

	/**
	 * Writes and read the alignment to and from XML.
	 * @author dmyersturnbull
	 * TODO Should use {@link AlignmentTools#toConciseAlignmentString(Map)}
	 */
	static class AlignmentFunctionAdapter extends XmlAdapter<String, Map<Integer, Integer>> {

		@Override
		public Map<Integer, Integer> unmarshal(String v) throws Exception {
			Map<Integer, Integer> map = new HashMap<Integer, Integer>();
			String[] comps = v.split(";");
			for (String comp : comps) {
				String[] parts = comp.split("=");
				map.put(Integer.parseInt(parts[0]), Integer.parseInt(parts[1]));
			}
			return map;
		}

		@Override
		public String marshal(Map<Integer, Integer> v) throws Exception {
			StringBuilder sb = new StringBuilder();
			for (Map.Entry<Integer, Integer> entry : v.entrySet()) {
				sb.append(entry.getKey() + "=" + entry.getValue() + ";");
			}
			return sb.toString();
		}
		
	}

	/**
	 * Constructs an AFPChain from this alignment mapping.
	 */
	public AFPChain buildAfpChain(Atom[] ca1, Atom[] ca2) throws StructureException {
		AFPChain afpChain = AlignmentTools.createAFPChain(ca1, ca2, new ResidueNumber[] {}, new ResidueNumber[] {});
		AlignmentTools.replaceOptAln(afpChain, ca1, ca2, getSimpleFunction());
		return afpChain;
	}
	
}
