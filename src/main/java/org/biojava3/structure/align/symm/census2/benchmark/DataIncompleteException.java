package org.biojava3.structure.align.symm.census2.benchmark;

public class DataIncompleteException extends RuntimeException {

	private static final long serialVersionUID = 7865520194833163469L;
	private String scopId;
	
	public String getScopId() {
		return scopId;
	}

	public DataIncompleteException(String scopId) {
		super();
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, String message) {
		super(message);
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, Throwable cause) {
		super(cause);
		this.scopId = scopId;
	}

	public DataIncompleteException(String scopId, String message, Throwable cause) {
		super(message, cause);
		this.scopId = scopId;
	}

	
}
