

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.util.ArrayList;

public class FullTreeTransform {

	private File dataFile;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		FullTreeTransform ftt = new FullTreeTransform();
		double[][] lambdaFullTree = ftt.getLambdaFullTree();
		double[][] llrFullTree = ftt.getLlrFullTree();
		long[] positions = ftt.getPositions();
		int[] snpBranches = ftt.getSNPbranches();
		double[] snpLambda = ftt.getSNPlambda();
		double[] snpLlr = ftt.getSNPllr();
		double[] snpR2 = ftt.getSNPr2();
		double[] snpCorrelations = ftt.getSNPcorrelations();
		int[] snpIndex = ftt.getSNPindex();
		String[] snpAncentry = ftt.getSNPancestry();
		String[] snpNames = ftt.getSNPnames();
		int[] treeMaxBranches = ftt.getTreeMaxBranches();
		double[] treeMaxLambda = ftt.getTreeMaxLambda();
		double[] treeMaxLlr = ftt.getTreeMaxLlr();
		double[] snpFreqs = ftt.getSnpFreqs();;
		
		return;
	}

	public FullTreeTransform(String[] args){

	}

	public FullTreeTransform(String x){

	}

	public FullTreeTransform(){

	}


	public double[][] getLlrFullTree(){

		this.dataFile = new File("llrFullTrees.dat");
		//use buffering
		ArrayList<Double[]> llrFullTrees = new ArrayList<Double[]>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			llrFullTrees = (ArrayList<Double[]>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int max = 0;

		for (int i = 0; i < llrFullTrees.size(); i++){
			if (llrFullTrees.get(i).length > max){
				max = llrFullTrees.get(i).length;
			}
		}
		double[][] llrResults = new double[llrFullTrees.size()][max];

		for (int i = 0; i < llrFullTrees.size(); i++){
			Double[] tempArray = new Double[max];
			if(llrFullTrees.get(i) != null){
				tempArray = llrFullTrees.get(i);
			}
			for (int j = 0; j < tempArray.length; j++){
				//				if(tempArray[i] != null){
				llrResults[i][j] = tempArray[j].doubleValue();
				//				}
			}
		}

		return llrResults;
	}

	public double[][] getLambdaFullTree(){

		this.dataFile = new File("lambdaFullTrees.dat");
		//use buffering
		ArrayList<Double[]> lambdaFullTrees = new ArrayList<Double[]>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			lambdaFullTrees = (ArrayList<Double[]>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int max = 0;

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			if (lambdaFullTrees.get(i).length > max){
				max = lambdaFullTrees.get(i).length;
			}
		}

		double[][] lambdaResults = new double[lambdaFullTrees.size()][max];

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			Double[] tempArray = new Double[max];
			if (lambdaFullTrees.get(i) != null){
				tempArray = lambdaFullTrees.get(i);
			} 
			for (int j = 0; j < tempArray.length; j++){
				if(tempArray[j] != null){
					lambdaResults[i][j] = tempArray[j].doubleValue();
				}
			}
		}

		return lambdaResults;
	}
	
	public double[][] getFreqsFullTree(){

		this.dataFile = new File("freqFullTrees.dat");
		//use buffering
		ArrayList<Double[]> lambdaFullTrees = new ArrayList<Double[]>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			lambdaFullTrees = (ArrayList<Double[]>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int max = 0;

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			if (lambdaFullTrees.get(i).length > max){
				max = lambdaFullTrees.get(i).length;
			}
		}

		double[][] lambdaResults = new double[lambdaFullTrees.size()][max];

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			Double[] tempArray = new Double[max];
			if (lambdaFullTrees.get(i) != null){
				tempArray = lambdaFullTrees.get(i);
			} 
			for (int j = 0; j < tempArray.length; j++){
				if(tempArray[j] != null){
					lambdaResults[i][j] = tempArray[j].doubleValue();
				}
			}
		}

		return lambdaResults;
	}

	public double[] getSNPllr(){
		this.dataFile = new File("snpLlr.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSnpFreqs(){
		this.dataFile = new File("snpFreqs.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSnpFreqsAtBranch(){
		this.dataFile = new File("snpFreqsAtBranch.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public double[] getSNPAncestralllr(){
		this.dataFile = new File("snpAncestralLlr.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public int[] getSNPbranches(){
		this.dataFile = new File("snpBranches.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}
	
	public int[] getSNPAncestralbranches(){
		this.dataFile = new File("snpAncestralBranches.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}
	
	public int[] getSNPindex(){
		this.dataFile = new File("snpIndex.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}

	public double[] getSNPlambda(){
		this.dataFile = new File("snpLambda.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPAncestrallambda(){
		this.dataFile = new File("snpAncestralLambda.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	

	public double[] getSNPr2(){
		this.dataFile = new File("snpR2.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPcorrelations(){
		this.dataFile = new File("snpCorrelations.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPAncestralcorrelations(){
		this.dataFile = new File("snpAncestralCorrelations.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public String[] getSNPancestry(){
		this.dataFile = new File("snpAncestry.dat");
		//use buffering
		ArrayList<String> snpLlrArray = new ArrayList<String>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<String>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		String[] result = new String[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i);
		}

		return result;
	}
	
	public String[] getSNPnames(){
		this.dataFile = new File("snpNames.dat");
		//use buffering
		ArrayList<String> snpLlrArray = new ArrayList<String>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<String>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		String[] result = new String[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i);
		}

		return result;
	}

	public double[] getTreeMaxLlr(){
		this.dataFile = new File("treeMaxLlr.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getTreeMaxFreqs(){
		this.dataFile = new File("treeMaxFreqs.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public double[] getTreeMaxLambda(){
		this.dataFile = new File("treeMaxLambda.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public int[] getTreeMaxBranches(){
		this.dataFile = new File("treeMaxBranch.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}

	public long[] getPositions(){
		this.dataFile = new File("positions.dat");
		//use buffering
		ArrayList<Long> snpBranchesArray = new ArrayList<Long>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Long>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		long[] result = new long[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).longValue();
		}

		return result;
	}
	
	//
	// Inverse!
	//
	
	public double[][] getLlrFullTreeI(){

		this.dataFile = new File("llrFullTreesI.dat");
		//use buffering
		ArrayList<Double[]> llrFullTrees = new ArrayList<Double[]>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			llrFullTrees = (ArrayList<Double[]>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int max = 0;

		for (int i = 0; i < llrFullTrees.size(); i++){
			if (llrFullTrees.get(i).length > max){
				max = llrFullTrees.get(i).length;
			}
		}
		double[][] llrResults = new double[llrFullTrees.size()][max];

		for (int i = 0; i < llrFullTrees.size(); i++){
			Double[] tempArray = new Double[max];
			if(llrFullTrees.get(i) != null){
				tempArray = llrFullTrees.get(i);
			}
			for (int j = 0; j < tempArray.length; j++){
				//				if(tempArray[i] != null){
				llrResults[i][j] = tempArray[j].doubleValue();
				//				}
			}
		}

		return llrResults;
	}

	public double[][] getLambdaFullTreeI(){

		this.dataFile = new File("lambdaFullTreesI.dat");
		//use buffering
		ArrayList<Double[]> lambdaFullTrees = new ArrayList<Double[]>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			lambdaFullTrees = (ArrayList<Double[]>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int max = 0;

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			if (lambdaFullTrees.get(i).length > max){
				max = lambdaFullTrees.get(i).length;
			}
		}

		double[][] lambdaResults = new double[lambdaFullTrees.size()][max];

		for (int i = 0; i < lambdaFullTrees.size(); i++){
			Double[] tempArray = new Double[max];
			if (lambdaFullTrees.get(i) != null){
				tempArray = lambdaFullTrees.get(i);
			} 
			for (int j = 0; j < tempArray.length; j++){
				if(tempArray[j] != null){
					lambdaResults[i][j] = tempArray[j].doubleValue();
				}
			}
		}

		return lambdaResults;
	}

	public double[] getSNPllrI(){
		this.dataFile = new File("snpLlrI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPAncestralllrI(){
		this.dataFile = new File("snpAncestralLlrI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public int[] getSNPbranchesI(){
		this.dataFile = new File("snpBranchesI.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}
	
	public int[] getSNPindexI(){
		this.dataFile = new File("snpIndexI.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}

	public double[] getSNPlambdaI(){
		this.dataFile = new File("snpLambdaI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPAncestrallambdaI(){
		this.dataFile = new File("snpAncestralLambdaI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public double[] getSNPr2I(){
		this.dataFile = new File("snpR2I.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPcorrelationsI(){
		this.dataFile = new File("snpCorrelationsI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public double[] getSNPAncestralcorrelationsI(){
		this.dataFile = new File("snpAncestralCorrelationsI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}
	
	public String[] getSNPancestryI(){
		this.dataFile = new File("snpAncestryI.dat");
		//use buffering
		ArrayList<String> snpLlrArray = new ArrayList<String>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<String>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		String[] result = new String[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i);
		}

		return result;
	}
	
	public String[] getSNPnamesI(){
		this.dataFile = new File("snpNamesI.dat");
		//use buffering
		ArrayList<String> snpLlrArray = new ArrayList<String>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<String>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		String[] result = new String[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i);
		}

		return result;
	}

	public double[] getTreeMaxLlrI(){
		this.dataFile = new File("treeMaxLlrI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public double[] getTreeMaxLambdaI(){
		this.dataFile = new File("treeMaxLambdaI.dat");
		//use buffering
		ArrayList<Double> snpLlrArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpLlrArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpLlrArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpLlrArray.get(i).doubleValue();
		}

		return result;
	}

	public int[] getTreeMaxBranchesI(){
		this.dataFile = new File("treeMaxBranchI.dat");
		//use buffering
		ArrayList<Integer> snpBranchesArray = new ArrayList<Integer>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Integer>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		int[] result = new int[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).intValue() + 1;
		}

		return result;
	}

	public double[] getPositionsI(){
		this.dataFile = new File("positionsI.dat");
		//use buffering
		ArrayList<Double> snpBranchesArray = new ArrayList<Double>();
		try{
			InputStream ifile = new FileInputStream(this.dataFile);
			InputStream buffer = new BufferedInputStream(ifile);
			ObjectInput input = new ObjectInputStream(buffer);

			snpBranchesArray = (ArrayList<Double>) input.readObject();

			input.close();
			buffer.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		double[] result = new double[snpBranchesArray.size()];

		for (int i = 0; i < result.length; i++){
			result[i] = snpBranchesArray.get(i).doubleValue();
		}

		return result;
	}
}
