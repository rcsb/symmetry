package org.biojava3.structure.codec;

import static org.rcsb.codec.CodecConstants.*;


import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Bond;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.rcsb.codec.BitEncoder;

/**
 * StructureToBinary writes/reads the "Atom" section of a PDB/mmCIF file to/from a binary files.
 * 
 * Format:
 * magic number "PDBb" (4 bytes)
 * version string (4 bytes)
 * 
 * The following records have the general format:
 * 
 * Record id, record length in number of bytes, data ..
 *    
 * Byte number  0   1   2   3   4   5   6   7   8   9
 *            +---+---+---+---+---+---+---+---+---+---+--
 *            |id | record length | data ..
 *            +---+---+---+---+---+---+---+---+---+---+--
 * 
 *            | N |             4 | model count (int = 4 bytes)
 *            
 *            | M |             4 | chains count in model (int = 4 bytes)
 *            
 *            | C |             8 | chain id (4 bytes), group count (int = 4 bytes)
 *            
 *            | G |            12 | short group record: group name (3 bytes), group number (int = 4 bytes), insertion code (1 byte), atom count (int = 4 bytes)
 *            
 *            | G |            24 | long group record: group name (3 bytes), group number (4 bytes), insertion code (1 byte), atom count (4 bytes), 
 *                                  xOffset(4 bytes), yOffset (4 bytes), zOffset (4 bytes)
 *            
 *            | A |            17 | short atom record: atom name with PDB spacing (4 bytes), element string (2 bytes), 
 *                                  alternative location indicator (1 byte), deltaX (short = 2 bytes), deltaY (short = 2 bytes), deltaZ (short = 2 bytes), 
 *                                  occupancy (short = 2 bytes), temperature factor (short = 2 bytes)
 *            
 *            | A |            27 | long atom record: atom name with PDB spacing (4 bytes), element string (2 bytes), 
 *                                  alternative location indicator (1 byte), x (float = 4 bytes), y (4 bytes), z (4 bytes), 
 *                                  occupancy (4 bytes), temperature factor (4 bytes)
 *                                  
 * for short atom record:  x = deltaX * 0.001 + xOffset; deltaY * 0.001 + yOffset; deltaZ * 0.001 + zOffset; 
 *                         o = occupancy * 0.001, b = temperature factor * 0.001;
 *            
 * Record order in file:
 * 
 *  magic number
 *  version
 *  N
 *  M model 1
 *  C chain 1
 *    G long group record 1    OR     G short group record 1
 *      A short atom record 1           A long atom record 1
 *      A short atom record 2           A long atom record 1
 *      ...
 *    G group 2
 *      ...
 *  C chain 2
 *  ..
 *  M model 2
 *  ..
 * 
 *
 * @author Peter Rose
 *
 */
public class StructureEncoderImpl1 extends StructureEncoder {
	private DataOutputStream outStream = null;
	
//	long byteCount = 0;
//	long totalCount = 0;
//	private int bCount = 0;
	private long time = 0;
	private int atomCount = 0;
	private int xyzSize = 0;
	private int bFactorSize = 0;
	private boolean useBfactor = false;
	private boolean useOccupancy = false;
	private boolean homogeneousModels = true;
	private List<Integer> groupArray = null;
	private List<int[]> bondInfo = null;
	private int[] intBuffer = new int[1];

	public StructureEncoderImpl1(DataOutputStream outStream) {
		this.outStream = outStream;
	}
	
	public void encode(Structure structure) throws IOException {
		atomCount = 0;
		groupArray = new ArrayList<Integer>();
		bondInfo =  new ArrayList<int[]>();
		writeInfo(structure);
		writeAtomInfo(structure);
		outStream.writeByte(END);
	}
	
	public void writeInfo(Structure structure) throws IOException {
		preCheck(structure);	
		outStream.writeByte(STRUCTURE);
		outStream.writeByte(5);
		outStream.writeInt(structure.nrModels());	
		outStream.writeBoolean(homogeneousModels);
		
		LinkedHashMap<String, Integer> groupMap = new LinkedHashMap<String, Integer>();
		LinkedHashMap<String, Integer> sequenceMap = new LinkedHashMap<String, Integer>();

		int modelCount = structure.nrModels();
		if (homogeneousModels) {
			modelCount = 1;
		}
			
		int groupIndex = 0;
		int sequenceIndex = 0;
		
		for (int m = 0 ; m < modelCount; m++) {
			List<Chain> chains = structure.getChains(m);
			outStream.writeByte(MODEL);
			outStream.writeByte(4);
			outStream.writeInt(chains.size());
			
			for (Chain c: chains) {
				int groupNumber = 0;	
				
				System.out.println(c.getChainID() + ": " + c.getSeqResSequence());
				Integer seqIndex = sequenceMap.get(c.getSeqResSequence());
				
				if (seqIndex == null) {
					sequenceMap.put(c.getSeqResSequence(), sequenceIndex);
					seqIndex = sequenceIndex;
					sequenceIndex++;
					outStream.writeByte(SEQUENCE);
					outStream.writeInt(c.getSeqResLength());
					writeFixedLengthString(c.getSeqResLength(), c.getSeqResSequence());
				}
				
				outStream.writeByte(CHAIN);
				outStream.writeByte(12);
				outStream.writeInt(seqIndex);
				writeFixedLengthString(4, c.getChainID());		
		//		List<Group> groups = c.getAtomGroups();
				List<Group> groups = c.getSeqResGroups();
				outStream.writeInt(groups.size());
				
				
				for (Group g: groups) {
					StringBuilder buffer = new StringBuilder();
					buffer.append(toFixedLengthString(3, g.getPDBName()));
					ResidueNumber rn = g.getResidueNumber();
					Character insCode = rn.getInsCode() == null ? ' ' : rn.getInsCode().charValue();
					buffer.append(insCode.toString());
					
					List<Atom> atoms = g.getAtoms();
					atomCount += atoms.size();
				
					for (Atom a: atoms) {
						buffer.append(a.getFullName());
						buffer.append(toFixedLengthString(2, a.getElement().name()));
						buffer.append(a.getAltLoc().toString());
					}
					
					String gInfo = buffer.toString();
					Integer gIndex = groupMap.get(gInfo);
				
					if (gIndex == null) {
						byte flags = 0;
						if (g instanceof AminoAcid) {
							flags |= AMINO_ACID;
							Atom head = null;
							try {
								head = ((AminoAcid) g).getN(); // outch, this has a nasty side effect of adding atoms
							} catch (StructureException e) {
							} 
							if (head != null) {
								flags |= HEAD;
							}
							Atom tail = null;
							try {
								tail = ((AminoAcid) g).getC();
							} catch (StructureException e) {
							} 
							if (tail != null) {
								flags |= TAIL;
							}	
						}
						if (g instanceof NucleotideImpl) {
							flags |= NUCLEOTIDE;
							Atom head = null;
							head = ((NucleotideImpl) g).getP();
							if (head != null) {
								flags |= HEAD;
							}
							Atom tail = null;
							tail = ((NucleotideImpl) g).getO3Prime();
							if (tail != null) {
								flags |= TAIL;
							}	
						}
						// this may return amino for an non-std amino acid
						if (g.getType().equals("HETATM")) {
							flags |= HETATOM;
						}
					
						gIndex = groupIndex;
						groupMap.put(gInfo, gIndex);
						groupIndex++;
						outStream.writeByte(GINFO);
						outStream.writeInt(gInfo.length()+ atoms.size()*4+3);
	//					System.out.println("write GINFO: " + (gInfo.length()+ atoms.size()*4+3));
						outStream.writeShort(atoms.size());
						outStream.writeByte(flags);
						outStream.write(gInfo.getBytes());
						int[] bondList = createBondList(g);
						for (int bl: bondList) {
						    outStream.writeShort(bl);
						}
						bondInfo.add(bondList);
					}
		
					groupArray.add(gIndex);			
					outStream.writeByte(GROUP);

					if (groupNumber == 0 || rn.getSeqNum() != groupNumber + 1) {
						// beginning of a chain or group numbers are discontinuous, save number
						groupNumber = rn.getSeqNum();
						outStream.writeByte(8);
						outStream.writeInt(gIndex);
						outStream.writeInt(groupNumber);
					} else {
						// sequential numbering, don't save value
						groupNumber++;
						outStream.writeByte(4);
						outStream.writeInt(gIndex);
					}
				}
			}
		}
		if (homogeneousModels) {
			atomCount *= structure.nrModels();
		}
//		System.out.println("groupArray: " + groupArray.size());
//		System.out.println("homogeneous: " + homogeneousModels);
//		System.out.println("useBfactor: " + useBfactor);
//		System.out.println("useOccupancy: " + useOccupancy);
//		System.out.println("models: " + structure.nrModels());
	}

	public void writeAtomInfo(Structure structure) throws IOException {
		List<Integer> xyz = new ArrayList<Integer>((int)((float)atomCount*3.6));
		int xyzType = 4;
		
		List<Integer> bFactor = null;
		List<Integer> occupancy = null;
		
		int bFactorType = 4;
		if (useBfactor) {
		   bFactor = new ArrayList<Integer>((int)((float)atomCount*1.2));
		}
		if (useOccupancy) {
		   occupancy = new ArrayList<Integer>((int)((float)atomCount));
		}
		
		int[] x = new int[1024];
		int[] y = new int[1024];
		int[] z = new int[1024];
		int[] b = new int[1024];
		
		int groupCount = 0;
		
		for (int m = 0 ; m < structure.nrModels(); m++) {
			if (homogeneousModels) {
			    groupCount = 0;
			}
			
			for (Chain c: structure.getChains(m)) {	
				// integer atom coordinates and b-factor for previous atom
				int xOffset = 0;
				int yOffset = 0;
				int zOffset = 0;
				int bOffset = 0;
				
				// integer atom coordinates and b-factor for previous 
				// link atom (Calpha for proteins, P for nucleic acids)
				int xLink = 0;
				int yLink = 0;
				int zLink = 0;
				int bLink = 0;
				boolean hasTail = false;
				
//				for (Group g: c.getAtomGroups()) {
					for (Group g: c.getSeqResGroups()) {
					int gIndex = groupArray.get(groupCount);
					int[] bondList = bondInfo.get(gIndex);
		
					List<Atom> atoms = g.getAtoms();
					System.out.println("atomCount: " + atoms.size());
					Atom tail = null;
					Atom head = null;
					boolean isAminoAcid = false;
					boolean isNucleotide = false;
					if (g instanceof AminoAcid) {
						head = null;
						try {
							head = ((AminoAcid) g).getN();
							
						} catch (StructureException e) {
						} 
						try {
							tail = ((AminoAcid) g).getC();
						} catch (StructureException e) {
						} 	
						isAminoAcid = true;
					} else if (g instanceof NucleotideImpl) {
						tail = ((NucleotideImpl) g).getO3Prime(); 
						head = ((NucleotideImpl) g).getP();
						isNucleotide = true;
					}
				
					if (atoms.size() > 1024) {
						x = new int[atoms.size()];
						y = new int[atoms.size()];
						z = new int[atoms.size()];
						b = new int[atoms.size()];
					}
					
					long t1 = System.nanoTime();
					for (int k = 0; k < atoms.size(); k++) {
						Atom a = atoms.get(k);
					
						x[k] = (int)Math.round(a.getX() * XYZ_PRECISION);
						y[k] = (int)Math.round(a.getY() * XYZ_PRECISION);
						z[k] = (int)Math.round(a.getZ() * XYZ_PRECISION);
						b[k] = (int)Math.round(a.getTempFactor() * B_PRECISION);
			
						int deltaX = 0;
						int deltaY = 0;
						int deltaZ = 0;
						int deltaB = 0;
						int distance = 0;
				
						if (bondList[k] >=0) {
							// use bond length from reference bond
							int reference = bondList[k];
							distance = bondList[atoms.size()+k];
							// calculate difference to previous bonded atom
							deltaX = x[k] - x[reference];
							deltaY = y[k] - y[reference];
							deltaZ = z[k] - z[reference];
							deltaB = b[k] - b[reference];
						} else if (k == 0 && head != null && hasTail) {
							// use standard polymer bond to previous residue, if it has a tail atom
							if (isAminoAcid) {
								distance = PEPTIDE_BOND_LENGTH;
							} else if (isNucleotide) {
								distance = NUCLEOTIDE_BOND_LENGTH;
							}
							// Calculate difference between link atom and link atom from previous group
							deltaX = x[k] - xLink;
							deltaY = y[k] - yLink;
							deltaZ = z[k] - zLink;
							deltaB = b[k] - bLink;
						} else {
							// calculate difference to previous atom
						
							deltaX = x[k] - xOffset;
							deltaY = y[k] - yOffset;
							deltaZ = z[k] - zOffset;
							deltaB = b[k] - bOffset;
						}
						
//						if (Math.abs(a.getZ() - 123.062) < 0.001) {
//							System.out.println("Atom: " + a);
//							System.out.println("k: " + k);
//							System.out.println("distance: " + distance);
//							System.out.println("Link atom: " + tail + " dist: " + distance);
//							System.out.println("offset: " + xOffset+","+yOffset+","+zOffset+","+bOffset);
//							System.out.println("delta: " + deltaX+","+deltaY+","+deltaZ+","+deltaB);
//							System.out.println("xyzb: " + x[k]+","+y[k]+","+z[k]+","+b[k]);
//							
//							
//							System.out.println("link: " + xLink+","+yLink+","+zLink);
//							System.out.println("------------------------------------");
//						}

						xyzType = encodeCoords(deltaX, deltaY, deltaZ, distance, xyzType, xyz);

						if (useBfactor) {
							bFactorType = addVarIntBfactor(deltaB, bFactorType, bFactor);
						}

						if (useOccupancy) {
							occupancy.add(Math.round((float)a.getOccupancy() * B_PRECISION));
						}			
					
						xOffset = x[k];
						yOffset = y[k];
						zOffset = z[k];
						bOffset = b[k];
					}
					time += System.nanoTime()-t1;
			
					if (tail != null) {
						hasTail = true;
						xLink = (int)Math.round(tail.getX() * XYZ_PRECISION);
						yLink = (int)Math.round(tail.getY() * XYZ_PRECISION);
						zLink = (int)Math.round(tail.getZ() * XYZ_PRECISION);
						bLink = (int)Math.round(tail.getTempFactor() * B_PRECISION);
					} else {
						hasTail = false;
						xLink = 0;
						yLink = 0;
						zLink = 0;
						bLink = 0;
					}	
					groupCount++;
				}

			}
		}
		
		bondInfo = null;
		groupArray = null;
		
		if (useBfactor) {
			outStream.writeByte(BFACTOR);
			outStream.writeInt(bFactorSize);
//			for (int bf: bFactor) {
//				System.out.println("b: " + bf);
//			}
			writeVarIntArray(bFactor);
		}
		if (useOccupancy) {
			outStream.writeByte(OCCUPANCY);
			outStream.writeInt(atomCount*2);
			for (int o: occupancy) {
				outStream.writeShort(o);
			}
		}
		
		outStream.writeByte(COORD);
//		System.out.println("writing coord record: " + xyzSize);
		outStream.writeInt(xyzSize);
		writeVarIntArray(xyz);
		
//		System.out.println("atomTime: " + (int)(time/1E6));
//		System.out.println("atomCount: " + atomCount);
//		System.out.println("xyzCount: " + xyz.size());
//		System.out.println("writeCoordlen: " + xyzSize);
//		System.out.println("bFactorSize: " + bFactorSize);
	}
	
	/**
	 * Encodes the x, y, z coordinates of an atom. When a bond distance is specified, it tries to encode x,y,z into a single 4 byte
	 * integer. Otherwise it encodes the x,y,z coordinates into 3 2-byte shorts or 3 4-byte integer, depending on the range of the x, y, z values.
	 * @param deltaX
	 * @param deltaY
	 * @param deltaZ
	 * @param distance
	 * @param intType
	 * @param array
	 * @return
	 * @throws IOException
	 */
	private int encodeCoords(int deltaX, int deltaY, int deltaZ, int distance, int intType, List<Integer> array) throws IOException {
		// if a standard bond distance is available, try encoding delta coordinates into a 4 byte integer
		if (distance != 0 && BitEncoder.toInt(distance, deltaX, deltaY, deltaZ, intBuffer)) {
			if (intType != 5) {
				array.add(getMarker(intType, 5));
				xyzSize += intType;
				intType = 5;
			}
			array.add(intBuffer[0]);
			xyzSize += 4;
		} else {	
			intType = addVarInt(deltaX, intType, array);
			intType = addVarInt(deltaY, intType, array);
			intType = addVarInt(deltaZ, intType, array);
		}
		return intType;
	}
	
	private int addVarInt(int value, int intType, List<Integer> array) throws IOException {
		int size = 2;	
		if (value < BYTE2_MIN_VALUE || value > BYTE2_MAX_VALUE) {
			size = 4;
		}
		if (value > BYTE4_MAX_VALUE) {
		   throw new IOException("addVarInt exceeds range");
		}

		switch (size) {
		case 2: 
			if (intType != 2) {
				array.add(getMarker(intType, 2));
				xyzSize += Math.min(intType, 4);
				intType = 2;
			}	
			break;
		case 4:
			if (intType != 4) {
				array.add(getMarker(intType, 4));
				xyzSize += Math.min(intType, 4);
				intType = 4;
			}
			break;
		}
		array.add(value);
		xyzSize += size;
		
		return intType;
	}
	
	private int getMarker(int prevType, int curType) {
		switch (prevType) {
		case 2: 
			return BYTE2_MAX_VALUE + curType;
		case 4:
			return BYTE4_MAX_VALUE + curType;
		case 5:
			return BYTE4_MAX_VALUE + curType;
		}
		return 8;
	}
	
	private void writeVarIntArray(List<Integer> values) throws IOException {
		int size = 4; // first element is always a 4-byte integer
		
		for (int i = 0; i < values.size(); i++) {
			int v = values.get(i);
			switch (size) {
			
			case 2:
				switch (v) {
				case BYTE2_MARKER4:
				case BYTE2_MARKER5: 
					outStream.writeShort(v);
					size = 4;
					v = values.get(++i);
					break;
				}
				break;
		
			case 4:
				switch (v) {
				case BYTE4_MARKER2: 
					outStream.writeInt(v);
					size = 2;
					v = values.get(++i);
					break;
				case BYTE4_MARKER5: 
					outStream.writeInt(v);
					v = values.get(++i);
					size = 4;
					break;
				}
			}
		
			if (size == 2) {
				outStream.writeShort(v);
			} else {
				outStream.writeInt(v);
			}
		}
	}

	
	private int addVarIntBfactor(int value, int intType, List<Integer> bFactors) {
		int size = 2;
		if (value < BYTE2_MIN_VALUE || value > BYTE2_MAX_VALUE) {
			size = 4;
		}

		switch (size) {
		case 2: 
			if (intType == 4) {
				bFactors.add(getMarker(intType, 2));
				bFactorSize += intType;
				intType = 2;
			}	
			break;
		
		case 4:
			if (intType == 2) {
				bFactors.add(getMarker(intType, 4));
				bFactorSize += intType;
				intType = 4;
			}
			break;
		}
		bFactors.add(value);
		bFactorSize += size;
		return intType;
	}

	private int[] createBondList(Group g) {
		List<Atom> atoms = g.getAtoms();
		int n = atoms.size();
		int[] bondList = new int[2*n];
		Arrays.fill(bondList, -1);

		for (int i = 0; i < n; i++) {
			Atom a = atoms.get(i);
			for (Bond b: a.getBonds()) {
				Atom other = b.getOther(a);
				int index = atoms.indexOf(other);
				if (index == -1) {
					continue;
				}
				if (index > i && bondList[index] == -1) {
					bondList[index] = i;
					try {
						bondList[n+index] = (int)Math.round(Calc.getDistance(a, other)* XYZ_PRECISION);
					} catch (StructureException e) {
					}
				}
			}
		}

		return bondList;
	}

	private String toFixedLengthString(int length, String string) throws IOException {
		if (string.length() == length) {
			return string;
		} else {
			StringBuilder sb = new StringBuilder(length);
			for (int i = 0, n = length - string.length(); i < n; i++) {
				sb.append(" ");
			}
			sb.append(string);
			return sb.toString();
		}
	}
	
	private void writeFixedLengthString(int length, String string) throws IOException {
		if (string.length() == length) {
			outStream.write(string.getBytes());
		} else {
			for (int i = 0, n = length - string.length(); i < n; i++) {
				outStream.write(BLANK);
			}
			outStream.write(string.getBytes());
		}
	}
	
	private void preCheck(Structure structure) {
		useBfactor = false;
		useOccupancy = false;
		homogeneousModels = true;
		
		// check first model for b-factor and occupancy
		for (Chain c: structure.getChains()) {	
			for (Group g: c.getAtomGroups()) {
				for (Atom a: g.getAtoms()) {
					if (a.getTempFactor() != 0.0) {
						useBfactor = true;
					}
					if (a.getOccupancy() != 1.0) {
						useOccupancy = true;
					}
				}
			}
		}
	
		// check if models are homogeneous
		for (int m = 1 ; m < structure.nrModels(); m++) {
			List<Chain> referenceChains = structure.getChains();
			List<Chain> chains = structure.getChains(m);
			
			if (referenceChains.size() != chains.size()) {
				homogeneousModels = false;
			}
			for (int i = 0; i < chains.size(); i++) {	
				Chain r = referenceChains.get(i);
				Chain c = chains.get(i);
				if (r.getAtomLength() != c.getAtomLength()) {
					homogeneousModels = false;
				}
				if (! r.getChainID().equals(c.getChainID())) {
					homogeneousModels = false;
				}
		
			    List<Group> referenceGroups = r.getAtomGroups();
			    List<Group> groups = c.getAtomGroups();
			    if (referenceGroups.size() != groups.size()) {
			    	homogeneousModels = false;
			    }
			    int groupCount = Math.min(referenceGroups.size(), groups.size());
			   
				for (int j = 0; j < groupCount; j++) {
					Group g = groups.get(j);
					Group rg = referenceGroups.get(j);

					if (homogeneousModels) {
						if (rg.size() != g.size()) {
							homogeneousModels = false;
						}
						if (! rg.getPDBName().equals(g.getPDBName())) {
							homogeneousModels = false;
						}
						if (rg.getResidueNumber().getSeqNum() != g.getResidueNumber().getSeqNum()) {
							homogeneousModels = false;
						}
					}

					List<Atom> referenceAtoms = rg.getAtoms();
					List<Atom> atoms = g.getAtoms();
					if (referenceAtoms.size() != atoms.size()) {
						homogeneousModels = false;
					}
					int atomCount = Math.min(referenceAtoms.size(), atoms.size());
					for (int k = 0; k < atomCount; k++) {
						Atom a = atoms.get(k);

						if (a.getTempFactor() != 0.0) {
							useBfactor = true;
						}
						if (a.getOccupancy() != 1.0) {
							useOccupancy = true;
						}
						if (homogeneousModels) {
							Atom ra = referenceAtoms.get(k);
							if (! a.getName().equals(ra.getName())) {
								homogeneousModels = false;
							}
						}
					}
				}
				// further checks necessary?
				if (! homogeneousModels && useBfactor && useOccupancy) {
					return;
				}
			}
		}
	}
	
}
