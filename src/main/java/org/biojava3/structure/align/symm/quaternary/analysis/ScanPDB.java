package org.biojava3.structure.align.symm.quaternary.analysis;


public class ScanPDB {

	public static void main(String args[]) {
		int threads = 1;
		for (int i = 0; i < threads; i++) {
			System.out.println("Thread: " + i);
			new Thread(new ScanPdbForQuarternarySymmetryNew(threads, i)).start();
		}
	}
}
