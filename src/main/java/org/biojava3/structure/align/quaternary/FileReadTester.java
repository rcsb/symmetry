package org.biojava3.structure.align.quaternary;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;

public class FileReadTester {

	public static void main(String[] args) throws InterruptedException {
		String path = "ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/divided/ft/1ft8.pdb.gz";
	//	String path = "ftp://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt";
	//	String path = "ftp://pdb.protein.osaka-u.ac.jp/pub/pdb/data/biounit/coordinates/divided/ft/1ft8.pdb.gz";
	// String path = "ftp://ftp.sun.com/pub/test.txt";
	//		String path = "ftp://pdb.protein.osaka-u.ac.jp/pub/pdb/data/biounit/coordinates/divided/u5/1u5b.pdb1.gz";
	//	String path = "ftp://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data/biounit/coordinates/divided/ft/1ft8.pdb.gz";
		for (int i = 0; i < 100; i++) {
			System.out.println("file exists: " + fileExists(path) + " " + i);
			Thread.sleep(500);
		}
	}

     public static boolean fileExists(String path) {
        URL url = null;
		try {
			url = new URL(path);
		} catch (MalformedURLException e) {
			e.printStackTrace();
		}
		InputStream is = null;
		URLConnection con = null;
    	try {
    		con = url.openConnection();
    //		con.setUseCaches(false);
		} catch (IOException e) {
			System.err.println("Got exception "+e);
			e.printStackTrace();
		}

		try {
			is=con.getInputStream();
		} catch (IOException e) {
			System.err.println("Got 2  exception "+e);
			e.printStackTrace();	
		}
		
		if (is == null) {
			return false;
		} else {
			try {
				is.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return true;
		}
     }
}
