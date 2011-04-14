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
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.changeux.IdentifyAllSymmetries;

public class ScanSCOPForSymmetry {
	public static void main(String[] args){
		ScanSCOPForSymmetry me = new ScanSCOPForSymmetry();

		me.scanSCOP("/Users/ap3/WORK/PDB/", true);
	}

	private void scanSCOP(String pdbFilePath, boolean isSplit) {

		AtomCache cache = new AtomCache(pdbFilePath, isSplit);
		ScopInstallation scop = new ScopInstallation(pdbFilePath);

		List<ScopDescription>superfamilies = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("got " + superfamilies.size() + " superfamilies");

		int fragmentLength = 8;
		
		int count = 0;
		int withSymm = 0;
		Map<Character,Integer> classStats = new HashMap<Character, Integer>();
		Map<Character,Integer> totalStats = new HashMap<Character, Integer>();
		for (ScopDescription superfamily : superfamilies){
			Character scopClass = superfamily.getClassificationId().charAt(0);
			
			if ( scopClass > 'g')
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
				StringBuffer str = new StringBuffer();
				if ( alternatives.size() > 0) {
					AFPChain afpChain = alternatives.get(0);
					if ( IdentifyAllSymmetries.isSignificant(afpChain)) {
						withSymm++;
						isSymmetric = true;
						
						
						
						     		
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
				} else {
					str.append("\t   ");  
					str.append("\t   ");  
					str.append("\t   ");  
				}
				
				System.out.print("#");
				if ( isSymmetric )
					System.out.print("* ");
				else
					System.out.print("  ");
				System.out.println(superfamily.getClassificationId() + "\t" + name1 + "\t" + alternatives.size() + "\t"+ count + "\t" + withSymm+ "\t" +
						String.format("%.2f",(withSymm/(float)count))  + "\t" + str.toString() +"\t" + superfamily.getDescription() + "\t" );

				trackStats(totalStats,scopClass,1);
				if ( isSymmetric) {
					trackStats(classStats,scopClass,1);
				}
				if ( withSymm > 5)
					break;
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
