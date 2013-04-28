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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;

import org.biojava.bio.structure.align.util.SynchronizedOutFile;

public class ResultConverter {

	public static String NEWLINE = System.getProperty("line.separator");

	/**
	 * Makes an HTML census file from an XML census file.
	 * @param args
	 * <ol>
	 * <li>The input XML file</li>
	 * <li>The name/path of the HTML file to create; cannot exist (just for safety)</li>
	 * </ol>
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: ResultConverter input-xml-file output-html-file");
			return;
		}
		File xmlIn = new File(args[0]);
		File htmlOut = new File(args[1]);
//		if (htmlOut.exists()) throw new IllegalArgumentException("File " + htmlOut + " already exists; cannot overwrite. Delete it and run again.");
		Results results = Results.fromXML(xmlIn);
		String html = results.toHTML();
		BufferedWriter bw = new BufferedWriter(new FileWriter(htmlOut));
		bw.write(html);
		bw.close();
	}

	public static String getHTMLFooter() {
		return "</tbody></table></body></html>" + NEWLINE;
	}

	public static String getHTMLHeader() {
		StringWriter file = new StringWriter();
		file.write("<!DOCTYPE html>" + NEWLINE);
		file.write("<html xmlns=\"http://www.w3.org/1999/xhtml\" lang=\"en\" xml:lang=\"en\">" + NEWLINE);
		file.write("\t<head>" + NEWLINE);
		file.write("\t\t<title>Census results</title>" + NEWLINE);
		file.write("\t</head>" + NEWLINE);
		file.write("\t<body>" + NEWLINE);
		file.write("\t<table id=\"census\">" + NEWLINE);
		file.write("\t\t<thead>" + NEWLINE);
		file.write("\t\t\t<tr>" + NEWLINE);
		file.write("\t\t\t<th>index</th>" + NEWLINE);
		file.write("\t\t\t<th>Symmetry</th>" + NEWLINE);
		file.write("\t\t\t<th>SCOP superfamily</th>" + NEWLINE);
		file.write("\t\t\t<th>SCOP name</th>" + NEWLINE);
		file.write("\t\t\t<th>TM-score</th>" + NEWLINE);
		file.write("\t\t\t<th>Z-score</th>" + NEWLINE);
		file.write("\t\t\t<th>RMSD</th>" + NEWLINE);
		file.write("\t\t\t<th>alignment score</th>" + NEWLINE);
		file.write("\t\t\t<th>% ID</th>" + NEWLINE);
		file.write("\t\t\t<th>% SIM </th>" + NEWLINE);
		file.write("\t\t\t<th>domain length</th>" + NEWLINE);
		file.write("\t\t\t<th>alignment length</th>" + NEWLINE);
		file.write("\t\t\t<th>angle</th>" + NEWLINE);
		file.write("\t\t\t<th>order</th>" + NEWLINE);
		file.write("\t\t\t<th>symmetry unit</th>" + NEWLINE);
		file.write("\t\t\t<th>SCOP description</th>" + NEWLINE);
		file.write("\t\t</tr>" + NEWLINE);
		file.write("\t</thead>" + NEWLINE);
		file.write("\t<tbody>" + NEWLINE);
		return file.toString();

	}

	public static void printHTMLHeader(SynchronizedOutFile file) throws IOException {
		String header = getHTMLHeader();
		file.write(header);
	}

	public static String toHTML(Result r) {

		int count = r.getRank();
		Boolean isSymm = r.getIsSignificant();
		String name = r.getScopId();
		StringBuilder str = new StringBuilder();
		str.append("\t\t<tr>" + NEWLINE);

		addCol(str, String.valueOf(count), isSymm);
		addCol(str, isSymm? "*" : "", isSymm);
		addCol(str, "<a href=\"http://scop.berkeley.edu/sunid=" + r.getSunId() + "\">" + r.getClassification() + "</a>", isSymm);
		addCol(str, "<a href=\"/jfatcatserver/showSymmetry.jsp?id=" + name + "\">" + name + "</a>", isSymm);

		if (r.getAlignment() != null) {
			addCol(str, String.format("%.2f", r.getAlignment().getTmScore()), isSymm);
			addCol(str, String.format("%.2f", r.getAlignment().getzScore()), isSymm);
			addCol(str, String.format("%.2f", r.getAlignment().getRmsd()), isSymm);
			addCol(str, String.format("%.2f", r.getAlignment().getAlignScore()), isSymm);
			addCol(str, String.format("%.2f%%", r.getAlignment().getIdentity()*100.0), isSymm);
			addCol(str, String.format("%.2f%%", r.getAlignment().getSimilarity()*100.0), isSymm);
			addCol(str, String.format("%d", r.getAlignment().getAlignLength() + r.getAlignment().getGapLength()), isSymm);
			addCol(str, String.format("%d", r.getAlignment().getAlignLength()), isSymm);
		} else {
			str.append("\t\t\t<td></td><td></td><td></td><td></td><td></td><td></td><td></td>" + NEWLINE);
		}
		
		if (r.getAxis() != null) {
			addCol(str, String.format("%.1f", r.getAxis().getTheta()), isSymm);
		} else {
			str.append("\t\t\t<td></td>");
		}
		
		if (r.getOrder() != null) {
			addCol(str, String.format("%d", r.getOrder()), isSymm);
		}
		
		if (r.getProtodomain() != null) {
			addCol(str, r.getProtodomain(), isSymm);
		} else {
			str.append("\t\t\t<td></td>" + NEWLINE);
		}

		if (r.getDescription() != null) {
			addCol(str, r.getDescription(), isSymm);
		} else {
			str.append("\t\t\t<td></td>" + NEWLINE);
		}

		str.append("\t\t</tr>" + NEWLINE);

		return str.toString();
	}

	private static void addCol(StringBuilder str, String text, boolean isSymmetric) {
		if (isSymmetric) {
			str.append("\t\t\t<td class=\"sig\">");
		} else {
			str.append("\t\t\t<td>");
		}
		str.append(text);
		str.append("</td>" + NEWLINE);
	}

	private ResultConverter() {

	}

}
