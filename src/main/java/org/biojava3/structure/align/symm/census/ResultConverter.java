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
 * Created on Sep 30, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.structure.align.symm.census;

import java.io.IOException;
import java.io.StringWriter;

import org.biojava.bio.structure.align.util.SynchronizedOutFile;


public class ResultConverter {

	public static String newline = System.getProperty("line.separator");

	private ResultConverter(){

	}

	public static void printHTMLHeader(SynchronizedOutFile file) throws IOException {
		String header = getHTMLHeader();
		file.write(header);
	}
	
	public static String getHTMLHeader(){
		StringWriter file = new StringWriter();
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
		file.write("<th>symmetry unit</th>");
		file.write("<th>SCOP description</th>");


		file.write("</tr>"+newline);
		file.write("<tbody>"+newline);
		return file.toString();

	}


	public static String toHTML(CensusResult r){

		int count = r.getRank();
		Boolean isSymmetric = r.getIsSignificant();
		String name = r.getName();


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
		str.append("<a href=\"http://scop.berkeley.edu/sunid="+r.getSunid()+"\">");
		//str.append("<a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid="+name+"\">");
		str.append(r.getClassificationId());
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
		str.append(String.format("%.2f",r.getzScore()));

		str.append("</td><td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(String.format("%.2f",r.getRmsd()));
		if ( isSymmetric)
			str.append("</b>");
		str.append("</td><td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(String.format("%.2f",r.getTmScore()));
		if ( isSymmetric)
			str.append("</b>");
		str.append("</td><td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(String.format("%.2f",r.getAlignScore()));
		if ( isSymmetric)
			str.append("</b>");
		str.append("</td>");


		addHtmlColumn(str,String.format("%.2f", r.getIdentity()), false);				
		addHtmlColumn(str,String.format("%.2f", r.getSimilarity()), false);									
		addHtmlColumn(str,String.format("%d",r.getAligLength()), isSymmetric);
		addHtmlColumn(str,String.format("%d",r.getAligLength()), isSymmetric);
		addHtmlColumn(str,String.format("%.1f",r.getAngle()), isSymmetric);			
		addHtmlColumn(str,String.format("%d",r.getOrder()), isSymmetric);
		addHtmlColumn(str, r.getProtoDomain(), isSymmetric);

		str.append("<td>");
		str.append(r.getDescription() );
		str.append("</td>");

		str.append("</tr>");

		return str.toString();
	}

	private  static void addHtmlColumn(StringBuffer str, String text, boolean isSymmetric) {
		str.append("<td>");
		if ( isSymmetric)
			str.append("<b>");
		str.append(text);
		if ( isSymmetric)
			str.append("</b>");

	}

}
