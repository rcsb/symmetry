package org.biojava3.structure.dbscan;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.changeux.IdentifyAllSymmetries;

//use script.ScanSCOPFrSymmetry instead
@Deprecated 
public class ScanSCOPForSymmetry {
	public static void main(String[] args){
		ScanSCOPForSymmetry me = new ScanSCOPForSymmetry();

		me.scanSCOP();
	}

	private void scanSCOP() {

		AtomCache cache = new AtomCache();
		ScopInstallation scop = new ScopInstallation();

		List<ScopDescription>superfamilies = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("got " + superfamilies.size() + " superfamilies");

		
		printHTMLHeader();
		int fragmentLength = 8;
		
		int count = 0;
		int withSymm = 0;
		Map<Character,Integer> classStats = new HashMap<Character, Integer>();
		Map<Character,Integer> totalStats = new HashMap<Character, Integer>();
		for (ScopDescription superfamily : superfamilies){
			Character scopClass = superfamily.getClassificationId().charAt(0);
			
			if ( scopClass > 'f')
				continue;
			
			count++;
			//System.err.println(superfamily);
			int sunid = superfamily.getSunID();
			List<ScopDomain> familyMembers = scop.getScopDomainsBySunid(sunid);
			ScopDomain first = familyMembers.get(0);
			try {
				String name1 = first.getScopId();
				String name2 = first.getScopId();
				
				IdentifyAllSymmetries identifyer = new IdentifyAllSymmetries();
				identifyer.setMaxNrAlternatives(1);
				identifyer.setDisplayJmol(false);
				
				List<AFPChain> alternatives = identifyer.indentifyAllSymmetries(name1, name2, cache, fragmentLength, null);
				boolean isSymmetric = false;
				
				AFPChain afpChain = null;

				if ( alternatives.size() > 0) {
					afpChain = alternatives.get(0);
					if ( IdentifyAllSymmetries.isSignificant(afpChain)) {
						withSymm++;
						isSymmetric = true;
					}
				}
						     		
				StringBuffer str = printTabbedResult(afpChain, isSymmetric, superfamily,name1, count);
				System.out.println(str.toString());
				trackStats(totalStats,scopClass,1);
				if ( isSymmetric) {
					trackStats(classStats,scopClass,1);
				}
				//if ( withSymm > 5)
				//	break;
			} catch (Exception e){
				e.printStackTrace();
			}

		}

		System.out.println("===");
		System.out.println("Overall symmetry: " + String.format("%.2f",(withSymm/(float)count)) + "%" );
		System.out.println("Statistics for SCOP classes:");
		for (Character scopClass: totalStats.keySet()){
			Integer total = totalStats.get(scopClass);
			Integer symm  = classStats.get(scopClass);
			System.out.println("Class: " + scopClass + " " + String.format("%.2f",(symm/(float)total))  + "%");
		}

	}
	

	private void printHTMLHeader() {
		System.out.println("<table>");
		System.out.println("<tr><td>#<\td><td>Symnmetry?<\td><td>SCOP superfamily<\td><td>SCOP name<\td>");
		System.out.println("<td>total % symmetry<\td><td>Z-score<\td><td>RMSD<\td><td>TM-score<\td><td>alignment score<\td>");
		System.out.println("<td>% ID<\td><td>% similarity<\td>");
		System.out.println("<td>SCOP description<\td><\tr>");
		
	}

	private StringBuffer printResult( AFPChain afpChain, boolean isSignificant, ScopDescription superfamily, String name, int count) {
		StringBuffer str = new StringBuffer();
		
		str.append("#");
		if ( isSignificant )
			str.append("* ");
		else
			str.append("  ");
		System.out.println(superfamily.getClassificationId() + "\t" + name + "\t" );
		
		
		if ( afpChain != null){
						str.append(String.format("%.2f",afpChain.getProbability()));		
						str.append("\t");
						str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));						
						str.append("\t");
						str.append(String.format("%.2f",afpChain.getTMScore()));
						str.append("\t");
						str.append(String.format("%.2f",afpChain.getAlignScore()));
						
					} else {
						str.append("\t   ");  
						str.append("\t   ");  
						str.append("\t   ");  
					}

		
		
		str.append(superfamily.getDescription() + "\t" );


		return str;
		
				}
				
	private StringBuffer printHTMLResult(AFPChain afpChain, boolean isSymmetric, ScopDescription superfamily, String name, int count){
		
		StringBuffer str = new StringBuffer();

		str.append("<tr>");
		str.append("<td>" + count+"<\td>");
				if ( isSymmetric )
			str.append("<td><b>*</b><\td>");
				else
			str.append("<td>&nbsp;<\td>");
		str.append("<td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(superfamily.getClassificationId());
		if ( isSymmetric)
			str.append("</b>");
		str.append("<\td><td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(
				"<a href=\"/jfatcatserver/analyzeSymmetry.jsp?scop1="+name + "\">"+name+"</a>" );
		if ( isSymmetric)
			str.append("</b>");		
		str.append("<\td>");

		
		str.append("<td>");


		if ( afpChain != null){		
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getProbability()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("<\td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("<\td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getTMScore()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("<\td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getAlignScore()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("<\td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f", afpChain.getIdentity()));			
			if ( isSymmetric)
				str.append("</b>");
			str.append("<\td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f", afpChain.getSimilarity()));			
			if ( isSymmetric)
				str.append("</b>");
			
		} else {
			str.append("<\td><td>   ");  
			str.append("<\td><td>   ");  
			str.append("<\td><td>   ");  
			str.append("<\td><td>   ");
			str.append("<\td><td>   ");
				}
		str.append("<\td>");

		str.append("<td>");
		str.append(superfamily.getDescription() );
		str.append("<\td>");
		
		str.append("<\tr>");
		
		return str;
		}

	private StringBuffer printTabbedResult(AFPChain afpChain, boolean isSymmetric, ScopDescription superfamily, String name, int count){
		
		StringBuffer str = new StringBuffer();

		str.append("");
		str.append("" + count+"\t");
		if ( isSymmetric )
			str.append("*\t");
		else
			str.append("\t");
		str.append("");
	
		str.append(superfamily.getClassificationId());
	
		str.append("\t");
		
		str.append(
				"/jfatcatserver/analyzeSymmetry.jsp?scop1="+name   );
		
		str.append("\t");
		
		
	


		if ( afpChain != null){		
		
			str.append(String.format("%.2f",afpChain.getProbability()));
	
			str.append("\t");
			
			str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
	
			str.append("\t");

			str.append(String.format("%.2f",afpChain.getTMScore()));
		
			str.append("\t");

			str.append(String.format("%.2f",afpChain.getAlignScore()));
	
			str.append("\t");

			str.append(String.format("%.2f", afpChain.getIdentity()));			
	
			str.append("\t");

			str.append(String.format("%.2f", afpChain.getSimilarity()));			

			
		} else {
			str.append("\t   ");  
			str.append("\t   ");  
			str.append("\t   ");  
			str.append("\t   ");
			str.append("\t   ");
		}
		str.append("\t");
		
		str.append("");
		str.append(superfamily.getDescription() );
		str.append("\t");
		
		
		
		return str;
	}

	
	
	private void trackStats(Map<Character, Integer> totalStats,
			Character scopClass, int i) {
	
		Integer number = totalStats.get(scopClass);
		if ( number == null) {
			number = 0;
			
		}
		
		number += i;
		totalStats.put(scopClass, number);
		
	}
}
