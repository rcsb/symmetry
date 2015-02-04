/*
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
 * Created on 2013-03-22
 *
 */
package org.biojava.nbio.structure;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;

import org.biojava.nbio.structure.scop.ScopInstallation;

/**
 * A simple {@link ScopDatabase} for mocking so that tests won't break just because SCOP changes.
 * This design (specifically, having this extend ScopInstallation and then adding lines to the cla, des, and hie files) really isn't a great choice,
 * but it was trivial to write.
 * @author dmyerstu
 */
public class MockScopDatabase extends ScopInstallation {

	public static final String claFileName = "mock.dir.cla.scop.";
	public static final String desFileName = "mock.dir.des.scop.";
	public static final String FILESPLIT;

	public static final String hieFileName = "mock.dir.hie.scop.";
	public static final String NEWLINE;
	String defaultDownloadUrl = "";

	String defaultScopVersion = "mock";
	static {
		NEWLINE = System.getProperty("line.separator");
		FILESPLIT = System.getProperty("file.separator");
	}

	/**
	 * Creates an empty {@link ScopDatabase} for mocking. Use the methods available to add lines.
	 */
	public MockScopDatabase() {
		super();
		new File(claFileName).delete();
		new File(desFileName).delete();
		new File(hieFileName).delete();
		setScopVersion(defaultScopVersion);
		setScopDownloadURL(defaultDownloadUrl);
	}

	public void close() {
		new File(claFileName).delete();
		new File(desFileName).delete();
		new File(hieFileName).delete();
	}
	
	/**
	 * Examples:
	 * 
	 * <pre>
	 * 46456	cl	a	-	All alpha proteins
	 * 164742	px	a.1.1.1	d2gkma_	2gkm A:
	 * 88965	sp	a.1.1.1	-	Mycobacterium tuberculosis, HbO [TaxId: 1773]
	 * </pre>
	 * 
	 * @param line
	 */
	public void addClaLine(String line) {
		write(claFileName, line + NEWLINE);
	}

	/**
	 * 
	 * @param lines
	 * @see #addClaLine(String)
	 */
	public void addClaLines(String[] lines) {
		StringBuilder sb = new StringBuilder();
		for (String line : lines) {
			sb.append(line + NEWLINE);
		}
		write(claFileName, sb.toString());
	}

	/**
	 * Examples:
	 * 
	 * <pre>
	 * 46456	cl	a	-	All alpha proteins
	 * 164742	px	a.1.1.1	d2gkma_	2gkm A:
	 * 88965	sp	a.1.1.1	-	Mycobacterium tuberculosis, HbO [TaxId: 1773]
	 * </pre>
	 * 
	 * @param line
	 */
	public void addDesLine(String line) {
		write(desFileName, line + NEWLINE);
	}

	/**
	 * 
	 * @param lines
	 * @see #addDesLine(String)
	 */
	public void addDesLines(String[] lines) {
		StringBuilder sb = new StringBuilder();
		for (String line : lines) {
			sb.append(line + NEWLINE);
		}
		write(desFileName, sb.toString());
	}

	/**
	 * Examples:
	 * 
	 * <pre>
	 * 46456	cl	a	-	All alpha proteins
	 * 164742	px	a.1.1.1	d2gkma_	2gkm A:
	 * 88965	sp	a.1.1.1	-	Mycobacterium tuberculosis, HbO [TaxId: 1773]
	 * </pre>
	 * 
	 * @param line
	 */
	public void addHieLine(String line) {
		write(hieFileName, line + NEWLINE);
	}

	/**
	 * 
	 * @param lines
	 * @see #addHieLine(String)
	 */
	public void addHieLines(String[] lines) {
		StringBuilder sb = new StringBuilder();
		for (String line : lines) {
			sb.append(line + NEWLINE);
		}
		write(hieFileName, sb.toString());
	}

	private void write(String file, String string) {
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(file));
			bw.write(string);
			bw.flush();
			bw.close();
		} catch (IOException e) {
			throw new IllegalStateException("Can't write to mock database file " + file, e);
		}
	}

	@Override
	protected void downloadClaFile() throws FileNotFoundException, IOException {
		
	}

	@Override
	protected void downloadDesFile() throws FileNotFoundException, IOException {
		
	}

	@Override
	protected void downloadHieFile() throws FileNotFoundException, IOException {
		
	}

	public void addClaLineWithHie(String string) {
		addClaLine(string);
		String[] parts = string.split("\\s")[4].split(",");
		// ex: cl=46456,cf=47768,sf=47819,fa=69044,dm=69045,sp=140642,px=129717
		String[] lines = new String[8];
		String cl = null, cf = null, sf = null, fa = null, dm = null, sp = null, px = null;
		for (String part : parts) {
			if (part.startsWith("cl=")) {
				cl = part.substring(3);
			} else if (part.startsWith("cf")) {
				cf = part.substring(3);
			} else if (part.startsWith("sf")) {
				sf = part.substring(3);
			} else if (part.startsWith("fa")) {
				fa = part.substring(3);
			} else if (part.startsWith("cf")) {
				dm = part.substring(3);
			} else if (part.startsWith("cf")) {
				sp = part.substring(3);
			} else if (part.startsWith("cf")) {
				px = part.substring(3);
			}
		}
		lines[0] = "0\t-\t" + cl;
		lines[1] = cl + "\t0\t" + cf;
		lines[2] = cf + "\t" + cl + "\t" + sf;
		lines[3] = sf + "\t" + cf + "\t" + fa;
		lines[4] = fa + "\t" + sf + "\t" + dm;
		lines[5] = dm + "\t" + fa + "\t" + sp;
		lines[6] = sp + "\t" + dm + "\t" + px;
		lines[7] = px + "\t" + sp + "\t-";
		addHieLines(lines);
	}

}
