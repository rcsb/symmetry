/**
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Sep 9, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package script;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.changeux.IdentifyAllSymmetries;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.utils.SymmetryTools;

// use org.biojava3.structure.align.symm.ScanSCOPForSymmetry instead
@Deprecated
public class ScanSCOPForSymmetry {

	public static String newline = System.getProperty("line.separator");

	public static StructureAlignment getCeSymm(){
		
		
		
		CeSymm ceSymm = new CeSymm();
		StructureAlignmentFactory.addAlgorithm(ceSymm);
		CeParameters params = (CeParameters) ceSymm.getParameters();
		if ( params == null) {
			params = new CeParameters();
			ceSymm.setParameters(params);
		}



		// here how to change the aa subst matrix, SDM is the default matrix
		String matrixName = "PRLA000101";
		SubstitutionMatrix<AminoAcidCompound> sdm = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName);			
		params.setSubstitutionMatrix(sdm);
		//SubstitutionMatrix<AminoAcidCompound> max = SubstitutionMatrixHelper.getBlosum85();
		//params.setSubstitutionMatrix(max);		

		// we over-weight sequence
		params.setSeqWeight(2.0);

		ceSymm.setParameters(params);

		return ceSymm;
	}

	public static void main(String[] args){

		String path =  "/Volumes/Macintosh HD2/PDB/";
		
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);
		

		try {
			ScanSCOPForSymmetry me = new ScanSCOPForSymmetry();

			me.run();

		} catch (Exception e){
			e.printStackTrace();
		}

	}

	private void run() {
		AtomCache cache = new AtomCache();
		ScopDatabase scop = ScopFactory.getSCOP();

		List<ScopDescription>superfamilies = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("got " + superfamilies.size() + " superfamilies");

		String f = cache.getPath() + File.separator + "scopCensus.html";

		try {
			SynchronizedOutFile outFile = new SynchronizedOutFile(new File(f));


			printHTMLHeader(outFile);

			//int fragmentLength = 8;

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

					if ( name1.equals("ds046__"))
						continue;
					
					Atom[] ca1 = cache.getAtoms(name1);
					Atom[] ca2 = cache.getAtoms(name2);

					boolean isSymmetric = false;

					StructureAlignment ceSymm = getCeSymm();
					AFPChain afpChain = ceSymm.align(ca1,ca2);

					double angle = -1;
					
					if ( afpChain != null) {

						angle = SymmetryTools.getAngle(afpChain,ca1,ca2);
						
						if ( IdentifyAllSymmetries.isSignificant(afpChain)) {
							
							if ( angle > 20) {
							
								withSymm++;
								isSymmetric = true;
							}
						}
					}
					
					//StringBuffer str = printTabbedResult(afpChain, isSymmetric, superfamily,name1, count);
					StringBuffer str = printHTMLResult(afpChain, isSymmetric, superfamily,name1, count, angle);
					//System.out.println(str.toString());
					outFile.write(str.toString() + newline);
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
			outFile.write("</tbody></table>");

			System.out.println("===");
			System.out.println("Overall symmetry: " + String.format("%.2f",(withSymm/(float)count)) + "%" );
			System.out.println("Statistics for SCOP classes:");
			for (Character scopClass: totalStats.keySet()){
				Integer total = totalStats.get(scopClass);
				Integer symm  = classStats.get(scopClass);
				System.out.println("Class: " + scopClass + " " + String.format("%.2f",(symm/(float)total))  + "%");
			}
			
			outFile.close();
		} catch (Exception e){
			e.printStackTrace();
		}


	}


	

	private void printHTMLHeader(SynchronizedOutFile file) throws IOException {
		file.write("<table id=\"census\">"+newline);
		file.write("<tr><th>index</th>");
		file.write("<th>Symmetry</th>");
		file.write("<th>SCOP superfamily</th>");
		file.write("<th>SCOP name</th>");
		file.write("<th>Z-score</th>");
		file.write("<th>RMSD</th>");
		file.write("<th>TM-score</th>");
		file.write("<th>alignment score</th>");
		file.write("<th>% ID</th>");
		file.write("<th>% SIM </th>");
		file.write("<th>domain length</th>");
		file.write("<th>alignment length</th>");
		file.write("<th>angle</th>");
		file.write("<th>order</th>");
		file.write("<th>SCOP description</th>");
		
		
		file.write("</tr>"+newline);
		file.write("<tbody>"+newline);

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

	private StringBuffer printHTMLResult(AFPChain afpChain, boolean isSymmetric, ScopDescription superfamily, String name, int count, double angle){

		StringBuffer str = new StringBuffer();

		str.append("<tr>");
		str.append("<td>" + count+"</td>");
		if ( isSymmetric )
			str.append("<td><b>*</b></td>");
		else
			str.append("<td>&nbsp;</td>");
		str.append("<td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append("<a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid="+name+"\">");
		str.append(superfamily.getClassificationId());
		str.append("</a>");
		if ( isSymmetric)
			str.append("</b>");
		str.append("</td><td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(
				"<a href=\"/jfatcatserver/showSymmetry.jsp?id="+name + "&matrix=sdm&seqWeight=2.0\">"+name+"</a>" );
		if ( isSymmetric)
			str.append("</b>");		
		str.append("</td>");


		str.append("<td>");


		if ( afpChain != null){		
			
			str.append(String.format("%.2f",afpChain.getProbability()));
			
			str.append("</td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("</td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getTMScore()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("</td><td>");
			if ( isSymmetric)
				str.append("<b>");
			str.append(String.format("%.2f",afpChain.getAlignScore()));
			if ( isSymmetric)
				str.append("</b>");
			str.append("</td>");
			
		
			addHtmlColumn(str,String.format("%.2f", afpChain.getIdentity()), false);				
			addHtmlColumn(str,String.format("%.2f", afpChain.getSimilarity()), false);									
			addHtmlColumn(str,String.format("%d",afpChain.getCa1Length()), isSymmetric);
			addHtmlColumn(str,String.format("%d",afpChain.getOptLength()), isSymmetric);
			addHtmlColumn(str,String.format("%.1f",angle), isSymmetric);
			
			int order = -1;
			try{
				order = CeSymm.getSymmetryOrder(afpChain);
			} catch (Exception e){
				e.printStackTrace();
			}
			addHtmlColumn(str,String.format("%d",order), isSymmetric);
			
			
			
		} else {
			str.append("</td><td>   ");  
			str.append("</td><td>   ");  
			str.append("</td><td>   ");  
			str.append("</td><td>   ");
			str.append("</td><td>   ");
			str.append("</td><td>   ");
			str.append("</td><td>   ");
			str.append("</td><td></td><td></td>");
		}
		

		str.append("<td>");
		str.append(superfamily.getDescription() );
		str.append("</td>");

		str.append("</tr>");

		return str;
	}

	private void addHtmlColumn(StringBuffer str, String text, boolean isSymmetric) {
		str.append("<td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(text);
		if ( isSymmetric)
			str.append("</b>");
		
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
				"/jfatcatserver/showSymmetry.jsp?id="+name   );

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

