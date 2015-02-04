/**
 * 
 */
package org.biojava.nbio.structure.codec;

import java.io.DataOutputStream;
import java.io.IOException;

import org.biojava.nbio.structure.Structure;

/**
 * @author Peter
 *
 */
public abstract class StructureEncoder {

	public static StructureEncoder getEncoder(int compressionMethod, DataOutputStream dataOutputStream) {
		if (compressionMethod == 1) {
			return new StructureEncoderImpl1(dataOutputStream);
		}
		return null;
	}
	
	public abstract void encode(Structure structure) throws IOException;
}
