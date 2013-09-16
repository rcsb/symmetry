package org.biojava3.structure.align.symm.census2.analysis;

import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

/**
 * A collection (map) of {@link StructureLigands}.
 * 
 * @author dmyerstu
 */
@XmlRootElement(name = "LigandList", namespace = "http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class LigandList implements Serializable {

	private static final long serialVersionUID = -6669870731175460575L;

	private static JAXBContext jaxbContext;
	private Map<String, StructureLigands> ligands;

	private String timestamp;

	static {
		try {
			jaxbContext = JAXBContext.newInstance(LigandList.class);
		} catch (Exception e) {
			throw new RuntimeException(e); // fatal
		}
	}

	public static LigandList fromXml(File file) throws IOException {

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();
			FileInputStream fis = new FileInputStream(file);
			LigandList ligands = (LigandList) un.unmarshal(fis);

			return ligands;

		} catch (JAXBException e) {
			throw new IOException(e);
		}

	}

	public LigandList() {
		ligands = new HashMap<String, StructureLigands>();
		timestamp = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date());
	}

	public boolean contains(String o) {
		return ligands.containsKey(o);
	}

	public Set<Entry<String, StructureLigands>> entrySet() {
		return ligands.entrySet();
	}

	public StructureLigands get(String key) {
		return ligands.get(key);
	}

	public Map<String, StructureLigands> getData() {
		return ligands;
	}

	public String getTimestamp() {
		return timestamp;
	}

	public boolean isEmpty() {
		return ligands.isEmpty();
	}

	public Set<String> keySet() {
		return ligands.keySet();
	}

	public void put(String key, StructureLigands value) {
		ligands.put(key, value);
	}

	public void add(String key, Ligand value) {
		if (!ligands.containsKey(key)) {
			ligands.put(key, new StructureLigands(key));
		}
		ligands.get(key).add(value);
	}

	public void putAll(Map<? extends String, ? extends Ligand> m) {
		for (Map.Entry<? extends String, ? extends Ligand> entry : m.entrySet()) {
			add(entry.getKey(), entry.getValue());
		}
	}

	public void setData(Map<String, StructureLigands> data) {
		this.ligands = data;
	}

	public void setTimestamp(String timestamp) {
		this.timestamp = timestamp;
	}

	public int size() {
		return ligands.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		return sb.toString();
	}

	public String toXml() throws IOException {

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		PrintStream ps = new PrintStream(baos);

		try {
			Marshaller m = jaxbContext.createMarshaller();
			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
			m.marshal(this, ps);
		} catch (JAXBException e) {
			throw new IOException(e);
		}

		return baos.toString();

	}

	public Collection<StructureLigands> values() {
		return ligands.values();
	}

	public void writeXmlToFile(File file) throws IOException {
		String xml = toXml();
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(file));
			bw.write(xml);
		} finally {
			if (bw != null)
				bw.close();
		}
	}

}
