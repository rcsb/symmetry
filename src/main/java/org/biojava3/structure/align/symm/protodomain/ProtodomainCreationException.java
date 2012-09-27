package org.biojava3.structure.align.symm.protodomain;

/**
 * An error in the creation of a Protodomain.
 * @author dmyersturnbull
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

	public ProtodomainCreationException(String protodomainString, String scopId, Throwable e, String reason) {
		super(makeMsg(protodomainString, scopId, reason), e);
		this.protodomainString = protodomainString;
		this.scopId = scopId;
	}

	public ProtodomainCreationException(String protodomainString, String scopId, Throwable e) {
		super(makeMsg(protodomainString, scopId, e.getMessage()), e);
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
