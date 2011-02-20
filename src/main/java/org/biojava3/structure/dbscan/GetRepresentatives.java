package org.biojava3.structure.dbscan;



import java.io.InputStream;

import java.net.URL;

import java.util.SortedSet;

import java.util.TreeSet;



import org.biojava.bio.structure.align.client.JFatCatClient;

import org.biojava.bio.structure.align.util.HTTPConnectionTools;



import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;

import org.rcsb.fatcat.server.PdbChainKey;



public class GetRepresentatives {

	static final String url = "http://www.rcsb.org/pdb/rest/representatives?cluster=40";



	public static SortedSet<PdbChainKey> getRepresentatives(){

		SortedSet<PdbChainKey> representatives = new TreeSet<PdbChainKey>();

		try {

			URL u = new URL(url);

			InputStream stream = HTTPConnectionTools.getInputStream(u,15000);



			String xml = null;



			if ( stream != null) {



				xml = JFatCatClient.convertStreamToString(stream);







				SortedSet<String> reps = RepresentativeXMLConverter.fromXML(xml);



				for (String s : reps){

					PdbChainKey k =  PdbChainKey.fromName(s);

					representatives.add(k);

				}

			}

		} catch(Exception e){

			e.printStackTrace();

		}

		return representatives;



	}

}


