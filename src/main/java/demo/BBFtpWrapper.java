package demo;

public class BBFtpWrapper {

	public static final String bbftpClient = "bbftpc";
	String cmdPath = null;
	public static void main(String[] args){
		String path = "/Users/ap3/Desktop/Downloads/bbftpPRO-7.0.4/";
		BBFtpWrapper bbftp = new BBFtpWrapper();
		
		bbftp.setCmdPath(path);
	}
	
	public void openConnection(String serverName){
		
	}
	
	public void closeConnection(){
		
	}
	
	/** Writes the content of a file to a newly created file at the target location
	 * 
	 * @param targetFileLocation
	 * @param fileContent
	 */
	public void writeToFile(String targetFileLocation, String fileContent){
		
	}
	
	
	public String getCmdPath() {
		return cmdPath;
	}
	public void setCmdPath(String cmdPath) {
		this.cmdPath = cmdPath;
	}
	
	
}
