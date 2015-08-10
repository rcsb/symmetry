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
 * Created on 2012-11-20
 *
 */

package org.biojava.nbio.structure.align.symm.protodomain;

/**
 * An error in the creation of a Protodomain.
 * 
 * @author dmyerstu
 * 
 */
public class ProtodomainCreationException extends Exception {

	private final String protodomainString;
	private final String scopId;

	private static String makeMsg(String protodomainString, String scopId, String reason) {
		return "Could not create protodomain " + protodomainString + " (scop ID " + scopId + "). Reason: " + reason;
	}

	public ProtodomainCreationException(String protodomainString, String scopId, String reason) {
		super(makeMsg(protodomainString, scopId, reason));
		this.protodomainString = protodomainString;
		this.scopId = scopId;
	}

	public ProtodomainCreationException(String protodomainString, String scopId, Throwable e) {
		super(makeMsg(protodomainString, scopId, e.getMessage()), e);
		this.protodomainString = protodomainString;
		this.scopId = scopId;
	}

	public ProtodomainCreationException(String protodomainString, String scopId, Throwable e, String reason) {
		super(makeMsg(protodomainString, scopId, reason), e);
		this.protodomainString = protodomainString;
		this.scopId = scopId;
	}

	public String getProtodomainString() {
		return protodomainString;
	}

	public String getScopId() {
		return scopId;
	}

}
