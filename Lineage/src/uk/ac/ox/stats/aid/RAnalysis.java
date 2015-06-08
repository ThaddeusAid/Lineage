package uk.ac.ox.stats.aid;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

public class RAnalysis {

	private File dataFile;
	private ArrayList<Integer> branches;
	private ArrayList<Double> llrArrayList;
	private ArrayList<Double[]> llrFullTreeArrayList;
	
	public void loadData (String dataFileString){
		this.dataFile = new File (dataFileString);
		
		try{
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
