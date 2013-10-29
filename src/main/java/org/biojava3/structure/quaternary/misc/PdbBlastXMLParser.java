/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.structure.quaternary.misc;


import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.biojava3.core.util.XMLHelper;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Peter Rose
 */
public class PdbBlastXMLParser {
    Document blastDoc = null;

    public PdbBlastXMLParser(InputStream blastXml) throws Exception {
        blastDoc = XMLHelper.inputStreamToDocument(blastXml);
    }

    public List<PdbBlastHit> parse(double maxEScore) throws Exception {
    	List<PdbBlastHit> hits = new ArrayList<PdbBlastHit>();

        ArrayList<Element> elementList = XMLHelper.selectElements(blastDoc.getDocumentElement(), "BlastOutput_iterations/Iteration[Iteration_hits]");

        for (Element element : elementList) {
            Element iterationHitsElement = XMLHelper.selectSingleElement(element, "Iteration_hits");
            ArrayList<Element> hitList = XMLHelper.selectElements(iterationHitsElement, "Hit");
            
            for (Element hitElement : hitList) {
                Element hitDefElement = XMLHelper.selectSingleElement(hitElement, "Hit_def");
                String hitDef = hitDefElement.getTextContent();
                Element hitLenElement = XMLHelper.selectSingleElement(hitElement, "Hit_len");
                String hitLen = hitLenElement.getTextContent();
                int hitLength = Integer.parseInt(hitLen);
               
                Element hithspsElement = XMLHelper.selectSingleElement(hitElement, "Hit_hsps");
                ArrayList<Element> hspList = XMLHelper.selectElements(hithspsElement, "Hsp");
                
                for (Element hspElement : hspList) {
                    Element evalueElement = XMLHelper.selectSingleElement(hspElement, "Hsp_evalue");
                    String value = evalueElement.getTextContent();
                    double evalue = Double.parseDouble(value);
                    if (evalue <= maxEScore) {
                    	Element lenElement = XMLHelper.selectSingleElement(hspElement, "Hsp_align-len");
                    	value = lenElement.getTextContent();
                    	int alignLen = Integer.parseInt(value);
                    	Element identityElement = XMLHelper.selectSingleElement(hspElement, "Hsp_identity");
                    	value = identityElement.getTextContent();
                    	int identity = Integer.parseInt(value);

                    	PdbBlastHit blastHit = parseBlastHit(hitDef);
                    	blastHit.setEvalue(evalue);
                    	blastHit.setIdentity((double)identity/alignLen);
                    	blastHit.setAlignmentLength(alignLen);
                    	blastHit.setHitLength(hitLength);

                    	hits.add(blastHit);
                    }
                }
            }
        }

        return hits;
    }

    /**
     * Parses content of Blast <Hit_def> element.
     * Example: 1IVP:1:A,B|pdbid|entity|chain(s)|sequence
     * @param hitString
     * @return PdbBlastHit
     */
	private PdbBlastHit parseBlastHit(String hitString) {
		int pos = hitString.indexOf("|");
		String[] items = hitString.substring(0, pos).split(":");
		
		PdbBlastHit hit = new PdbBlastHit();
		hit.setPdbId(items[0].trim());
		hit.setChainIds(Arrays.asList(items[2].split("[ ,]+")));
        
		return hit;
	}
	
	public static void sortBlastHits(List<PdbBlastHit> blastHits) {
		Collections.sort(blastHits, new Comparator<PdbBlastHit>() {
			@Override
			public int compare(PdbBlastHit o1, PdbBlastHit o2) {
				double ratio1 =  o1.getIdentity() * (1.0 - Math.abs(o1.getAlignmentLength() - o1.getHitLength())/(double)o1.getAlignmentLength());
				double ratio2 =  o2.getIdentity() *(1.0 - Math.abs(o2.getAlignmentLength() - o2.getHitLength())/(double)o2.getAlignmentLength());
                double val = ratio2 - ratio1;
//				System.out.println("Compare: " + val);
//				System.out.println(o1);
//		    	System.out.println(o2);
                return (int)Math.signum(val);
			}
		});
	}
}
