package org.biojava3.structure.align.symm.census2;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.annotation.adapters.XmlAdapter;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;

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
	
}
