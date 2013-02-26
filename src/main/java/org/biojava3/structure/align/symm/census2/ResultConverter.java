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
package org.biojava3.structure.align.symm.census2;

import java.io.IOException;
import java.io.StringWriter;

import org.biojava.bio.structure.align.util.SynchronizedOutFile;

public class ResultConverter {

	public static String newline = System.getProperty("line.separator");

	public static String getHTMLFooter() {
		return "</tbody></table></body></html>" + newline;
	}

	public static String getHTMLHeader() {
		StringWriter file = new StringWriter();
		file.write("<!DOCTYPE html>" + newline);
		file.write("<html xmlns=\"http://www.w3.org/1999/xhtml\" lang=\"en\" xml:lang=\"en\">" + newline);
		file.write("<head><title>Census results</title</head>" + newline);
		file.write("<body>" + newline);
		file.write("<table id=\"census\">" + newline);
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

		file.write("</tr>" + newline);
		file.write("<tbody>" + newline);
		return file.toString();

	}

	public static void printHTMLHeader(SynchronizedOutFile file) throws IOException {
		String header = getHTMLHeader();
		file.write(header);
	}

	public static String toHTML(Result r) {

		int count = r.getRank();
		Boolean isSymmetric = r.getIsSignificant();
		String name = r.getScopId();

		StringBuffer str = new StringBuffer();

		str.append("<tr>");
		str.append("<td>" + count + "</td>");
		if (isSymmetric) str.append("<td><strong>*</strong></td>");
		else str.append("<td>&nbsp;</td>");
		str.append("<td>");
		if (isSymmetric) str.append("<strong>");
		str.append("<a href=\"http://scop.berkeley.edu/sunid=" + r.getSunId() + "\">");
		str.append(r.getClassification());
		str.append("</a>");
		if (isSymmetric) str.append("</strong>");
		str.append("</td><td>");
		if (isSymmetric) str.append("<strong>");
		str.append("<a href=\"/jfatcatserver/showSymmetry.jsp?id=" + name + "&matrix=sdm&seqWeight=2.0\">" + name
				+ "</a>");
		if (isSymmetric) str.append("</strong>");
		str.append("</td>");
		str.append("<td>");
		str.append(String.format("%.2f", r.getAlignment().getzScore()));

		str.append("</td><td>");
		if (isSymmetric) str.append("<strong>");
		str.append(String.format("%.2f", r.getAlignment().getRmsd()));
		if (isSymmetric) str.append("</strong>");
		str.append("</td><td>");
		if (isSymmetric) str.append("<strong>");
		str.append(String.format("%.2f", r.getAlignment().getTmScore()));
		if (isSymmetric) str.append("</strong>");
		str.append("</td><td>");
		if (isSymmetric) str.append("<strong>");
		str.append(String.format("%.2f", r.getAlignment().getAlignScore()));
		if (isSymmetric) str.append("</strong>");
		str.append("</td>");

		addHtmlColumn(str, String.format("%.2f", r.getAlignment().getIdentity()), false);
		addHtmlColumn(str, String.format("%.2f", r.getAlignment().getSimilarity()), false);
		addHtmlColumn(str, String.format("%d", r.getAlignment().getAlignLength()), isSymmetric);
		addHtmlColumn(str, String.format("%d", r.getAlignment().getAlignLength()), isSymmetric);
		addHtmlColumn(str, String.format("%.1f", r.getAxis().getTheta()), isSymmetric);
		addHtmlColumn(str, String.format("%d", r.getOrder()), isSymmetric);
		addHtmlColumn(str, r.getProtodomain(), isSymmetric);

		str.append("<td>");
		str.append(r.getDescription());
		str.append("</td>");

		str.append("</tr>");

		return str.toString();
	}

	private static void addHtmlColumn(StringBuffer str, String text, boolean isSymmetric) {
		str.append("<td>");
		if (isSymmetric) str.append("<strong>");
		str.append(text);
		if (isSymmetric) str.append("</strong>");

	}

	private ResultConverter() {

	}

}
