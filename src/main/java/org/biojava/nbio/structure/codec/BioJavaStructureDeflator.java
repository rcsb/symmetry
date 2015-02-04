/**
 * 
 */
package org.biojava.nbio.structure.codec;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Structure;
import static org.rcsb.codec.CodecConstants.*;

/**
 * @author Peter
 *
 */
public class BioJavaStructureDeflator {
	private static final byte MAJOR_VERSION = 0;
	private static final byte MINOR_VERSION = 0;

	private int compressionMethod = 1;	
	long byteCount = 0;
	long totalCount = 0;
	
	private long fileSize = 0;
	private long fileSizeCompressed = 0;
	private long writeTime = 0;
	
	private DataOutputStream outStream = null;
	private File file = null;
	

	/**
	 * @return the fileSize
	 */
	public long getFileSize() {
		return fileSize;
	}

	/**
	 * @return the fileSizeCompressed
	 */
	public long getFileSizeCompressed() {
		return fileSizeCompressed;
	}

	/**
	 * @return the writeTime
	 */
	public long getWriteTime() {
		return writeTime;
	}

	public void deflate(Structure structure, String fileName, int compressionMethod) throws IOException {
		this.compressionMethod = compressionMethod;
		
		fileSize = 0;
		fileSizeCompressed = 0;
		writeTime = 0;
		
		long start = System.nanoTime();

		open(fileName);
		
		writeHeader();
		writeAtoms(structure);
		outStream.flush();
		
		close();
		fileSize = outStream.size();
		outStream.close();
		
		File file = new File(fileName);
		fileSizeCompressed = file.length();
		
		writeTime = System.nanoTime() - start;

	//	System.out.println("Wrote bytes: raw: " + fileSize +  " file size: " + fileSizeCompressed + " in "+ (writeTime/1E6) + " ms" );
	}
	
	public void append(Structure structure) throws IOException {			
		long start = System.nanoTime();
		
		writeHeader();
		writeAtoms(structure);
		outStream.flush();
		
		writeTime += start - System.nanoTime();

		fileSize = outStream.size();
	
	//	System.out.println("Appended bytes: raw: " + fileSize +  " file size: " + fileSizeCompressed + " in "+ (writeTime/1E6) + " ms" );
	}
	
	public void close() throws IOException {
		if (outStream != null) {
		   outStream.close();
		}
	}
	
	private void writeHeader() throws IOException {
		outStream.write(MAGIC_NUMBER.getBytes());
		outStream.writeByte(MAJOR_VERSION);
		outStream.writeByte(MINOR_VERSION);
		outStream.writeByte(compressionMethod);
	}
	
	private void open(String fileName) throws IOException {
		file = new File(fileName);
		file.createNewFile();
		file.setWritable(true);
	
		outStream = new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(file))));
	}
	
	private void writeAtoms(Structure structure) throws IOException {
		StructureEncoder encoder = StructureEncoder.getEncoder(compressionMethod, outStream);
		encoder.encode(structure);
	}
	
}
