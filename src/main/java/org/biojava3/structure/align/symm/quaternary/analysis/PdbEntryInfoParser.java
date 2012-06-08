package org.biojava3.structure.align.symm.quaternary.analysis;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava3.structure.align.symm.quaternary.analysis.PdbEntryInfo.ExperimentalMethod;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public final class PdbEntryInfoParser {
	private static String URL = "http://www.rcsb.org/pdb/rest/getBioAssemblies";
	private static Map<String, ExperimentalMethod> expMethodMap = new HashMap<String, ExperimentalMethod>();
	
	private static List<PdbEntryInfo> pdbEntryInfo = new ArrayList<PdbEntryInfo>();

	static {	
		expMethodMap.put("electronCrystallography", ExperimentalMethod.ELECTRON_CRYSTALLOGRAPHY);
		expMethodMap.put("electronMicroscopy", ExperimentalMethod.ELECTRON_MICROSCOPY);
		expMethodMap.put("solidStateNmr", ExperimentalMethod.FIBER_DIFFRACTION);
		expMethodMap.put("neutronDiffraction", ExperimentalMethod.NEUTRON_DIFFRACTION);
		expMethodMap.put("solidStateNmr", ExperimentalMethod.SOLID_STATE_NMR);
		expMethodMap.put("nmr", ExperimentalMethod.SOLUTION_NMR);
		expMethodMap.put("solutionScattering", ExperimentalMethod.SOLUTION_SCATTERING);
		expMethodMap.put("theoreticalModel", ExperimentalMethod.THEORETICAL_MODEL);
		expMethodMap.put("xray", ExperimentalMethod.X_RAY);
		
		parseXML();
	}
	
	private PdbEntryInfoParser() {}; // this class is a singleton, it cannot be instantiated

	public static List<PdbEntryInfo> getPdbEntryInfo() {
		return pdbEntryInfo;
	}

	private static final void parseXML() {
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-mm-dd");

		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = null;
		try {
			db = factory.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		}
		Document doc = null;
		try {
			doc = db.parse(URL);
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		doc.getDocumentElement().normalize();


		NodeList nodes = doc.getElementsByTagName("PDB");

		for(int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
				
			PdbEntryInfo info = new PdbEntryInfo();
			
			// parse experiment method
			List<String> methods = getAttributes(node, "method", "name");
			List<ExperimentalMethod> expMethods = Collections.synchronizedList(new ArrayList<ExperimentalMethod>(methods.size()));
			for (String exp: methods) {
				expMethods.add(expMethodMap.get(exp));
			}
			info.setExperimentalMethods(expMethods);
			
			// parse entity type
			List<String> entityTypes = getAttributes(node, "Entry", "type");	
			info.setEntityTypes(entityTypes);
			
			// parse pdbId
			NamedNodeMap map = node.getAttributes();
			String pdbId = null;
			Node nd = map.getNamedItem("structureId");
			if (nd == null) {
				continue;
			}
			pdbId = nd.getTextContent();
			info.setPdbId(pdbId);

			// parse bioassemblies
			int bioAssemblyCount = 0;
			nd = map.getNamedItem("bioAssemblies");
			if (nd != null) {
				bioAssemblyCount = Integer.parseInt(nd.getTextContent());
			}
			info.setBioAssemblyCount(bioAssemblyCount);

			// parse release date
			Date releaseDate = null;		
			nd = map.getNamedItem("release_date");
			if (nd != null) {
				try {
					releaseDate =  formatter.parse(nd.getTextContent());
				} catch (ParseException e) {
					e.printStackTrace();
				}
			}
			info.setReleaseDate(releaseDate);

			// parse resolution
			float resolution = Float.NaN;		
			nd = map.getNamedItem("resolution");
			if (nd != null) {
				resolution =  Float.parseFloat(nd.getTextContent());
			}
			info.setResolution(resolution);

			info.setPdbId(pdbId);
			pdbEntryInfo.add(info);
		}

	}

	private static List<String> getAttributes(Node node, String tagName, String attributeName) {
		ArrayList<String> list = new ArrayList<String>();
		
		if (node.getNodeType() == Node.ELEMENT_NODE) {
			Element e = (Element) node;
			NodeList methodsNodes = e.getElementsByTagName(tagName);
			if (methodsNodes != null) {
				for (int j = 0; j < methodsNodes.getLength(); j++) {
					Node methodsNode = methodsNodes.item(j);
					NamedNodeMap map = methodsNode.getAttributes();
					if (map != null) {
						Node nn = map.getNamedItem(attributeName);
						list.add(nn.getTextContent());
					}
				}
			}
		}
		
		list.trimToSize();
		return list;
	}
}
