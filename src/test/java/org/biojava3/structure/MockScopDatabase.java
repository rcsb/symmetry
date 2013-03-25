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
package org.biojava3.structure;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;

import org.biojava.bio.structure.scop.ScopInstallation;

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
		String remoteFilename = claFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getClaFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	@Override
	protected void downloadDesFile() throws FileNotFoundException, IOException {
		String remoteFilename = desFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getDesFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	@Override
	protected void downloadHieFile() throws FileNotFoundException, IOException {
		String remoteFilename = hieFileName + scopVersion + ".txt";
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getHieFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

}
