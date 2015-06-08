package uk.ac.ox.stats.aid;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.Semaphore;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class TreeAnalysis extends Thread {

	public Semaphore semaphore;

	private String[] currentTreeStringArray;
	private double[][] snps;
	private String[] referenceAllele;
	private String[] alternateAllele;
	private String[] ancestorAllele;
	private String[] snpNames;
	private int[] members;
	private boolean usingSMC;

	private long[] legendPositions;

	private int numberOfTaxa;
	private long range;
	private long currentPosition;
	private long endPosition;
	private double r2Threshold;
	private double llr;
	private long[] positions;
	private ArrayList<Integer> snpBranch;
	private ArrayList<Integer> snpIndex;
	private ArrayList<String> snpAncestor;
	private ArrayList<String> snpName;

	private Double[] llrResults;
	private Double[] lambdaResults;

	private ArrayList<Double> maxLambda;
	private ArrayList<Integer> maxRow;
	private ArrayList<Double> snpR2;

	private ArrayList<Integer> treeBranch;

	private ArrayList<Double> snpLlr;
	private ArrayList<Double> snpAncestralLlr;
	private ArrayList<Double> snpLambda;
	private ArrayList<Double> snpAncestralLambda;
	private ArrayList<Double> snpCorrelation;
	private ArrayList<Double> snpAncestralCorrelation;
	private ArrayList<Integer> snpAncestorBranch;

	private ArrayList<Double> treeLlr;
	private ArrayList<Double> treeLambda;

	private ArrayList<Long> snpPositions;

	private ArrayList<Double> snpFreqs;
	private ArrayList<Double> snpFreqsAtBranch;
	private Double[] branchFreqs;
	private ArrayList<Double> treeMaxFreqs;

	public TreeAnalysis(double[][] SNPs, String[] treeStringArray, long[] positions, int taxa, long RANGE, long currentPOSITION, double r2, Semaphore sem, String[] refAllele, String[] altAllele, String[] ancestorAllele, String[] snpN){
		this.snps = SNPs;
		this.currentTreeStringArray = treeStringArray;
		this.legendPositions = positions;
		this.numberOfTaxa = taxa;
		this.range = RANGE;
		this.currentPosition = currentPOSITION;
		this.r2Threshold = r2;
		this.semaphore = sem;

		this.referenceAllele = refAllele;
		this.alternateAllele = altAllele;
		this.ancestorAllele = ancestorAllele;
		this.snpNames = snpN;
		this.positions = positions;
		this.members = null;
		this.usingSMC = false;
	}

	public TreeAnalysis(double[][] SNPs, String[] treeStringArray, long[] positions, int taxa, long RANGE, long currentPOSITION, double r2, Semaphore sem, String[] refAllele, String[] altAllele, String[] ancestorAllele, String[] snpN, int[] memberTranslation){
		this.snps = SNPs;
		this.currentTreeStringArray = treeStringArray;
		this.legendPositions = positions;
		this.numberOfTaxa = taxa;
		this.range = RANGE;
		this.currentPosition = currentPOSITION;
		this.r2Threshold = r2;
		this.semaphore = sem;

		this.referenceAllele = refAllele;
		this.alternateAllele = altAllele;
		this.ancestorAllele = ancestorAllele;
		this.snpNames = snpN;
		this.positions = positions;
		this.members = memberTranslation;
		this.usingSMC = true;
	}

	@Override
	public synchronized void start() {
		super.start();
	}

	@Override
	public void run() {
		super.run();

		int[][] currentTreeMatrix = new int[0][0];
		int popSize = 0;

		// init arraylists
		snpBranch = new ArrayList<Integer>();
		snpIndex = new ArrayList<Integer>();
		snpAncestor = new ArrayList<String>();;
		snpName = new ArrayList<String>();

		maxLambda = new ArrayList<Double>();
		maxRow = new ArrayList<Integer>();
		snpR2 = new ArrayList<Double>();

		treeBranch = new ArrayList<Integer>();

		snpLlr = new ArrayList<Double>();
		snpLambda = new ArrayList<Double>();
		snpCorrelation = new ArrayList<Double>();

		this.snpAncestralLlr = new ArrayList<Double>();
		this.snpAncestralLambda = new ArrayList<Double>();
		this.snpAncestralCorrelation = new ArrayList<Double>();
		this.snpAncestorBranch = new ArrayList<Integer>();

		treeLlr = new ArrayList<Double>();
		treeLambda = new ArrayList<Double>();

		snpPositions = new ArrayList<Long>();

		this.snpFreqs = new ArrayList<Double>();
		this.snpFreqsAtBranch = new ArrayList<Double>();
		this.branchFreqs = new Double[0];
		this.treeMaxFreqs = new ArrayList<Double>();

		try{
			if(this.usingSMC == true){
				// if smc file
				currentTreeMatrix = this.makeNewickTree(currentTreeStringArray[3]);
				popSize = currentTreeMatrix.length + 1;
				this.endPosition = Long.parseLong(currentTreeStringArray[2]);
			} else {
				// if genecluster tree file
				currentTreeMatrix = this.makeTree(currentTreeStringArray[2].split(","));
				currentTreeStringArray[1] = currentTreeStringArray[1].trim();

				popSize = Integer.parseInt(currentTreeStringArray[1]);

				this.endPosition = this.currentPosition;
			}

			//			System.out.println("Start Position: " + this.currentPosition + " End Position: " + this.endPosition);
			this.numberOfTaxa = popSize; // not sure why I have both of these

			// look for evidence of selection
			AidZhantikinMyers zam = new AidZhantikinMyers(currentTreeMatrix, popSize, popSize - 1, 20, 0.2);
			zam.run();

			// try and associate LLR results to individual SNPs
			// yVectors[numberOfBranches][numberOfTaxa]
			double[][] yVectors = new double[this.numberOfTaxa - 1][this.numberOfTaxa];

			for (int branch = 0; branch < this.numberOfTaxa - 1; branch++){

				Arrays.fill(yVectors[branch], 0.0);
				ArrayList<Integer> currentMembers = this.getMembers(currentTreeMatrix, branch);

				if(this.members == null){
					for (int j = 0; j < currentMembers.size(); j++){
						yVectors[branch][currentMembers.get(j)] = 1.0;
					}
				} else {
					for (int j = 0; j < currentMembers.size(); j++){
						yVectors[branch][this.members[currentMembers.get(j)]] = 1.0;
					}
				}
			}

			double lowerLimit = this.currentPosition;
			double upperLimit = this.endPosition;
			int rangeStart = 0;
			int rangeEnd = 0;
			int currentSNPBranch = -1;

			// find the start of the SNP range
			if (this.legendPositions == null){
				this.legendPositions = new long[0];
			}
			for (int j = 0; j < this.legendPositions.length; j++){
				if (lowerLimit >= this.legendPositions[j]){
					rangeStart = j;
				} 

				if (upperLimit >= this.legendPositions[j]){
					rangeEnd = j;
				}
			}

			// match the SNPs to the branches
			PearsonsCorrelation pc = new PearsonsCorrelation();
			int[] branches = new int[rangeEnd - rangeStart + 1];
			Arrays.fill(branches,  -1);
			double[] correlations = new double[rangeEnd - rangeStart + 1];
			Arrays.fill(correlations, Double.NEGATIVE_INFINITY);

			int[] branchesI = new int[rangeEnd - rangeStart + 1];
			Arrays.fill(branchesI,  -1);
			double[] correlationsI = new double[rangeEnd - rangeStart + 1];
			Arrays.fill(correlationsI, Double.NEGATIVE_INFINITY);

			// wait for the threads to complete
			zam.join();

			this.llrResults = zam.getLlrResults();
			this.lambdaResults = zam.getLambdaResults();
			this.branchFreqs = zam.getBranchFreqs();

			for (int snps = rangeStart; snps <= rangeEnd; snps++){
				if (this.snps == null){
					break;
				}

				double[] snpVector = new double[this.numberOfTaxa];

				int numberOfAltAlleles = 0;

				// get the snp vector
				for ( int j = 0; j < this.numberOfTaxa; j++){
					snpVector[j] = this.snps[j][snps];
					numberOfAltAlleles += this.snps[j][snps];
				}

				String ancestry;

				// 1000 genomes phase 1 ancestry and phase 3 ? entries
				if(this.referenceAllele[snps].trim().equalsIgnoreCase(this.ancestorAllele[snps].trim())){
					ancestry = "a";
				} else if (this.alternateAllele[snps].trim().equalsIgnoreCase(this.ancestorAllele[snps].trim())){
					ancestry = "r";
				} else {
					ancestry = "u";
				}

				// 1000 genomes phase 3 ancestry
				if(this.ancestorAllele[snps].contains("|")){
					if(this.ancestorAllele[snps].contains("unknown")){
						ancestry = "u";
					} else if (this.ancestorAllele[snps].contains("insertion")){
						if(this.referenceAllele[snps].length() < this.alternateAllele[snps].length()){
							ancestry = "a";
						} else {
							ancestry = "r";
						}
					} else if(this.ancestorAllele[snps].contains("deletion")){
						if(this.referenceAllele[snps].length() > this.alternateAllele[snps].length()){
							ancestry = "a";
						} else {
							ancestry = "r";
						}
					} else {
						String[] currentArray = this.ancestorAllele[snps].split("|");
						
						if(this.referenceAllele[snps].trim().equalsIgnoreCase(currentArray[0].trim())){
							ancestry = "a";
						} else if (this.alternateAllele[snps].trim().equalsIgnoreCase(currentArray[0].trim())){
							ancestry = "r";
						} else {
							ancestry = "u";
						}
					}
				}
				
				// flip reference derived alleles so that 1 is derived
				if(ancestry.equals("r")){
					for(int i = 0; i < snpVector.length; i++){
						if(snpVector[i] == 0){
							snpVector[i] = 1;
						} else {
							snpVector[i] = 0;
						}
					}
				}

				// if site is monomorphic or singleton skip to next SNP
				if(numberOfAltAlleles == 1 || numberOfAltAlleles == 0 || numberOfAltAlleles == this.numberOfTaxa || numberOfAltAlleles == this.numberOfTaxa - 1){
					continue;
				}

				double[] currentCorrelations = new double[this.numberOfTaxa - 1];
				double[] currentR2 = new double[this.numberOfTaxa - 1];

				//System.out.println(ancestry);
				// correlate the current snp against all branches
				for (int j = 0; j < yVectors.length; j++){
					currentCorrelations[j] = pc.correlation(snpVector, yVectors[j]);
					currentR2[j] = (double)(float)Math.pow(currentCorrelations[j], 2.0);
				}

				double max = Double.NEGATIVE_INFINITY;
				double maxCorr = 0.0;
				int maxIndex = -1;
				double currentLlr = Double.NEGATIVE_INFINITY;

				double maxa = Double.NEGATIVE_INFINITY;
				double maxCorra = 0.0;
				int maxIndexa = -1;
				double currentLlra = Double.NEGATIVE_INFINITY;

				for (int j = 0; j < currentCorrelations.length; j++){ 

					if (ancestry.equals("u") && currentR2[j] > max){
						max = currentR2[j];
						maxIndex = j;
						maxCorr = currentCorrelations[j];

						currentLlr = this.llrResults[j];
					} else if (ancestry.equals("u") && currentR2[j] == max && currentLlr < this.llrResults[j]){
						max = currentR2[j];
						maxIndex = j;
						maxCorr = currentCorrelations[j];

						currentLlr = this.llrResults[j];
					} else {
						if (currentCorrelations[j] > 0 && currentR2[j] > max){
							max = currentR2[j];
							maxIndex = j;
							maxCorr = currentCorrelations[j];

							currentLlr = this.llrResults[j];

						}

						if (currentCorrelations[j] < 0 && currentR2[j] > maxa){
							maxa = currentR2[j];
							maxIndexa = j;
							maxCorra = currentCorrelations[j];

							currentLlra = this.llrResults[j];

						}
					}
				}


				branches[snps - rangeStart] = maxIndex;
				correlations[snps - rangeStart] = max;

				double freq = 0.0;

				for(int i = 0; i < snpVector.length; i++){
					freq += (double) snpVector[i];
				}

				freq /= (double) snpVector.length;

				this.snpFreqs.add(freq);

				if(maxIndex > -1){
					this.snpFreqsAtBranch.add(this.branchFreqs[maxIndex]);
				} else {
					this.snpFreqsAtBranch.add(-1.0);
				}
					int snpFreq = 0;

					for (int i = 0; i < snpVector.length; i++){
						snpFreq += snpVector[i];
					}
				
				int snpFreqI = snpVector.length - snpFreq;

				this.snpIndex.add(snps);
				currentSNPBranch = maxIndex;
				int currentSNPBrancha = maxIndexa;
				this.snpR2.add(max);
				this.snpCorrelation.add(maxCorr);
				this.snpAncestralCorrelation.add(maxCorra);
				this.snpName.add(this.snpNames[snps]);
				this.snpPositions.add(positions[snps]);

				if (ancestry.equals("r")){
					this.snpAncestor.add("Alt");
				} else if (ancestry.equals("a")){
					this.snpAncestor.add("Ref");
				} else {
					this.snpAncestor.add("Ukn");
				}

				if (snpFreq == 0 || snpFreq == snpVector.length){
					currentSNPBranch = this.numberOfTaxa - 1;
					this.snpR2.add(1.0);
					if(snpVector[0] == 1){
						this.snpCorrelation.add(-1.0);
					}
					if(snpVector[0] == 0){
						this.snpCorrelation.add(1.0);
					}
				}
				if(currentSNPBranch > -1 && currentSNPBranch < this.numberOfTaxa - 1 && currentSNPBrancha > -1 && currentSNPBrancha < this.numberOfTaxa - 1){
					this.snpLlr.add(llrResults[currentSNPBranch].doubleValue());
					this.snpLambda.add(this.lambdaResults[currentSNPBranch].doubleValue());

					this.snpAncestralLlr.add(llrResults[currentSNPBrancha].doubleValue());
					this.snpAncestralLambda.add(this.lambdaResults[currentSNPBrancha].doubleValue());

				} else {
					this.snpLlr.add(0.0);
					this.snpLambda.add(0.0);

					this.snpAncestralLlr.add(0.0);
					this.snpAncestralLambda.add(0.0);
				}

				this.snpBranch.add(currentSNPBranch);
				this.snpAncestorBranch.add(currentSNPBrancha);

				int col = -1;
				max = Double.NEGATIVE_INFINITY;

				for (int i = 0; i < this.llrResults.length; i++){
					if (max < this.llrResults[i].doubleValue()){
						col = i;
						max = this.llrResults[i].doubleValue();
					}
				}


				this.maxLambda.add(this.lambdaResults[col]);
				this.treeLambda.add(this.lambdaResults[col]);
				this.maxRow.add(col);
				this.treeLlr.add(max);
				this.treeMaxFreqs.add(this.branchFreqs[col]);
			}
		

		if (this.snps == null){
			this.snpIndex.add(null);
			this.snpR2.add(null);
			this.snpCorrelation.add(null);
			this.snpName.add(null);
			this.snpPositions.add(this.currentPosition);


			this.snpAncestor.add(null);

			this.snpBranch.add(null);

			int col = -1;
			double max = Double.NEGATIVE_INFINITY;

			for (int i = 0; i < this.llrResults.length; i++){
				if (max < this.llrResults[i].doubleValue()){
					col = i;
					max = this.llrResults[i].doubleValue();
				}
			}


			this.maxLambda.add(this.lambdaResults[col]);
			this.treeLambda.add(this.lambdaResults[col]);
			this.maxRow.add(col);
			this.treeLlr.add(max);
		}

		zam = null;

	} catch (Exception e){
		e.printStackTrace();
	}
	this.semaphore.release();
}


public int[][] makeTree(String[] tempStringArray) {

	if (this.numberOfTaxa == 0){
		this.numberOfTaxa = tempStringArray.length + 1;
	}

	int[][] currentTree = new int[this.numberOfTaxa - 1][3];

	for (int j = 0; j < tempStringArray.length; j++) {
		String[] branches = tempStringArray[j].split(" ");
		currentTree[j][0] = Integer.parseInt(branches[1]) - 1;
		currentTree[j][1] = Integer.parseInt(branches[2]) - 1;
		currentTree[j][2] = this.numberOfTaxa + j;
		//      System.out.println(this.treesArray[i][j][0] + "," + this.treesArray[i][j][1] + " " + this.treesArray[i][j][2]);
	}

	return currentTree;
}

public ArrayList<Integer> getMembers(int[][] tree, int branch) {
	ArrayList<Integer> members = new ArrayList<Integer>();

	if (tree[branch][0] < this.numberOfTaxa) {
		Integer member1 = tree[branch][0];
		members.add(member1);
	} else {
		members.addAll(this.getMembers(tree, tree[branch][0] - this.numberOfTaxa));
	}

	if (tree[branch][1] < this.numberOfTaxa) {
		Integer member2 = tree[branch][1];
		members.add(member2);
	} else {
		members.addAll(this.getMembers(tree, tree[branch][1] - this.numberOfTaxa));
	}

	return members;
}

public ArrayList<Integer> getBranch() {
	return snpBranch;
}

public double getLlr() {
	return llr;
}

public Double[] getLlrResults() {
	return llrResults;
}

public ArrayList<Double> getMaxLambda() {
	return maxLambda;
}

public ArrayList<Integer> getMaxRow() {
	return maxRow;
}

public Double[] getLambdaResults() {
	return lambdaResults;
}

public ArrayList<Integer> getMaxBranch() {
	return treeBranch;
}

public ArrayList<Integer> getTreeBranch() {
	return treeBranch;
}

public ArrayList<Double> getTreeLambda() {
	return treeLambda;
}

public ArrayList<Double> getTreeLlr() {
	return treeLlr;
}

public ArrayList<Double> getR2() {
	return snpR2;
}

public double getR2Threshold() {
	return r2Threshold;
}

public long getRange() {
	return range;
}

public ArrayList<Integer> getSnpBranch() {
	return snpBranch;
}

public ArrayList<Double> getSnpLambda() {
	return snpLambda;
}

public ArrayList<Double> getSnpLlr() {
	return snpLlr;
}

public double[][] getSnps() {
	return snps;
}

public ArrayList<Double> getSnpCorrelation() {
	return snpCorrelation;
}

public ArrayList<Integer> getSnpIndex(){
	return snpIndex;
}

public ArrayList<String> getSnpAncestor() {
	return snpAncestor;
}

public ArrayList<String> getSnpName() {
	return snpName;
}

//	public ArrayList<Double> getTreeLlrI() {
	//		return treeLlrI;
//	}

@Override
public ClassLoader getContextClassLoader() {
	// TODO Auto-generated method stub
	return super.getContextClassLoader();
}

@Override
public long getId() {
	// TODO Auto-generated method stub
	return super.getId();
}

@Override
public StackTraceElement[] getStackTrace() {
	// TODO Auto-generated method stub
	return super.getStackTrace();
}

public String[] getAlternateAllele() {
	return alternateAllele;
}

public String[] getAncestorAllele() {
	return ancestorAllele;
}

public double getCurrentPosition() {
	return currentPosition;
}

@Override
public State getState() {
	// TODO Auto-generated method stub
	return super.getState();
}

@Override
public UncaughtExceptionHandler getUncaughtExceptionHandler() {
	// TODO Auto-generated method stub
	return super.getUncaughtExceptionHandler();
}

public String[] getCurrentTreeStringArray() {
	return currentTreeStringArray;
}

//	public Double[] getLambdaResultsI() {
//		return lambdaResultsI;
//	}

public long[] getLegendPositions() {
	return legendPositions;
}

//	public Double[] getLlrResultsI() {
//		return llrResultsI;
//	}

//	public ArrayList<Double> getMaxLambdaI() {
//		return maxLambdaI;
//	}

//	public ArrayList<Integer> getMaxRowI() {
//		return maxRowI;
//	}

public int getNumberOfTaxa() {
	return numberOfTaxa;
}

public String[] getReferenceAllele() {
	return referenceAllele;
}

public Semaphore getSemaphore() {
	return semaphore;
}

//	public ArrayList<Double> getSnpCorrelationI() {
//		return snpCorrelationI;
//	}

//	public ArrayList<Double> getSnpLambdaI() {
//		return snpLambdaI;
//	}

//	public ArrayList<Double> getSnpLlrI() {
//		return snpLlrI;
//	}

public String[] getSnpNames() {
	return snpNames;
}

public ArrayList<Double> getSnpR2() {
	return snpR2;
}

//	public ArrayList<Double> getSnpR2I() {
//		return snpR2I;
//	}

//	public ArrayList<Integer> getTreeBranchI() {
//		return treeBranchI;
//	}

//	public ArrayList<Double> getTreeLambdaI() {
//		return treeLambdaI;
//	}

public int[][] makeNewickTree(String input){
	//		input = "((((((27:1002.805899[&&NHX:age=0.000000],(86:0.000000[&&NHX:age=0.000000],13:0.000000[&&NHX:age=0.000000])172:1002.805899[&&NHX:age=0.000000])107:0.000000[&&NHX:age=1002.805899],1:1002.805899[&&NHX:age=0.000000])101:0.000000[&&NHX:age=1002.805899],((63:0.000000[&&NHX:age=0.000000],32:0.000000[&&NHX:age=0.000000])126:0.000000[&&NHX:age=0.000000],(((94:0.000000[&&NHX:age=0.000000],((90:0.000000[&&NHX:age=0.000000],67:0.000000[&&NHX:age=0.000000])180:0.000000[&&NHX:age=0.000000],(51:0.000000[&&NHX:age=0.000000],40:0.000000[&&NHX:age=0.000000])102:0.000000[&&NHX:age=0.000000])134:0.000000[&&NHX:age=0.000000])188:0.000000[&&NHX:age=0.000000],(43:0.000000[&&NHX:age=0.000000],((((99:0.000000[&&NHX:age=0.000000],70:0.000000[&&NHX:age=0.000000])198:0.000000[&&NHX:age=0.000000],37:0.000000[&&NHX:age=0.000000])140:0.000000[&&NHX:age=0.000000],(85:0.000000[&&NHX:age=0.000000],55:0.000000[&&NHX:age=0.000000])110:0.000000[&&NHX:age=0.000000])170:0.000000[&&NHX:age=0.000000],36:0.000000[&&NHX:age=0.000000])147:0.000000[&&NHX:age=0.000000])171:0.000000[&&NHX:age=0.000000])159:0.000000[&&NHX:age=0.000000],((60:0.000000[&&NHX:age=0.000000],(78:0.000000[&&NHX:age=0.000000],(96:0.000000[&&NHX:age=0.000000],(30:0.000000[&&NHX:age=0.000000],21:0.000000[&&NHX:age=0.000000])165:0.000000[&&NHX:age=0.000000])192:0.000000[&&NHX:age=0.000000])156:0.000000[&&NHX:age=0.000000])120:0.000000[&&NHX:age=0.000000],7:0.000000[&&NHX:age=0.000000])169:0.000000[&&NHX:age=0.000000])143:0.000000[&&NHX:age=0.000000])127:1002.805899[&&NHX:age=0.000000])158:4361.040322[&&NHX:age=1002.805899],11:5363.846221[&&NHX:age=0.000000])119:2687.856147[&&NHX:age=5363.846221],((52:1002.805899[&&NHX:age=0.000000],(54:639.178343[&&NHX:age=0.000000],92:639.178343[&&NHX:age=0.000000])108:363.627557[&&NHX:age=639.178343])133:2559.449441[&&NHX:age=1002.805899],(22:1002.805899[&&NHX:age=0.000000],(19:1002.805899[&&NHX:age=0.000000],33:1002.805899[&&NHX:age=0.000000])149:0.000000[&&NHX:age=1002.805899])173:2559.449441[&&NHX:age=1002.805899])160:4489.447028[&&NHX:age=3562.255340])129:18918.895955[&&NHX:age=8051.702368],((75:1545.314509[&&NHX:age=0.000000],24:1545.314509[&&NHX:age=0.000000])184:16499.310954[&&NHX:age=1545.314509],(74:12061.808515[&&NHX:age=0.000000],(34:2354.701987[&&NHX:age=0.000000],((64:1545.314509[&&NHX:age=0.000000],((81:0.000000[&&NHX:age=0.000000],97:0.000000[&&NHX:age=0.000000])194:0.000000[&&NHX:age=0.000000],79:0.000000[&&NHX:age=0.000000])162:1545.314509[&&NHX:age=0.000000])105:809.387478[&&NHX:age=1545.314509],(88:395.449492[&&NHX:age=0.000000],((((39:0.000000[&&NHX:age=0.000000],8:0.000000[&&NHX:age=0.000000])155:0.000000[&&NHX:age=0.000000],4:0.000000[&&NHX:age=0.000000])121:0.000000[&&NHX:age=0.000000],((89:0.000000[&&NHX:age=0.000000],84:0.000000[&&NHX:age=0.000000])178:0.000000[&&NHX:age=0.000000],2:0.000000[&&NHX:age=0.000000])168:0.000000[&&NHX:age=0.000000])113:0.000000[&&NHX:age=0.000000],(((20:0.000000[&&NHX:age=0.000000],((41:0.000000[&&NHX:age=0.000000],(71:0.000000[&&NHX:age=0.000000],(58:0.000000[&&NHX:age=0.000000],((80:0.000000[&&NHX:age=0.000000],(28:0.000000[&&NHX:age=0.000000],45:0.000000[&&NHX:age=0.000000])111:0.000000[&&NHX:age=0.000000])150:0.000000[&&NHX:age=0.000000],23:0.000000[&&NHX:age=0.000000])179:0.000000[&&NHX:age=0.000000])116:0.000000[&&NHX:age=0.000000])142:0.000000[&&NHX:age=0.000000])163:0.000000[&&NHX:age=0.000000],6:0.000000[&&NHX:age=0.000000])181:0.000000[&&NHX:age=0.000000])157:0.000000[&&NHX:age=0.000000],((49:0.000000[&&NHX:age=0.000000],44:0.000000[&&NHX:age=0.000000])195:0.000000[&&NHX:age=0.000000],((68:0.000000[&&NHX:age=0.000000],(66:0.000000[&&NHX:age=0.000000],50:0.000000[&&NHX:age=0.000000])132:0.000000[&&NHX:age=0.000000])136:0.000000[&&NHX:age=0.000000],(((17:0.000000[&&NHX:age=0.000000],47:0.000000[&&NHX:age=0.000000])104:0.000000[&&NHX:age=0.000000],12:0.000000[&&NHX:age=0.000000])187:0.000000[&&NHX:age=0.000000],((31:0.000000[&&NHX:age=0.000000],(25:0.000000[&&NHX:age=0.000000],16:0.000000[&&NHX:age=0.000000])197:0.000000[&&NHX:age=0.000000])106:0.000000[&&NHX:age=0.000000],(38:0.000000[&&NHX:age=0.000000],((82:0.000000[&&NHX:age=0.000000],((72:0.000000[&&NHX:age=0.000000],83:0.000000[&&NHX:age=0.000000])182:0.000000[&&NHX:age=0.000000],(59:0.000000[&&NHX:age=0.000000],93:0.000000[&&NHX:age=0.000000])118:0.000000[&&NHX:age=0.000000])166:0.000000[&&NHX:age=0.000000])164:0.000000[&&NHX:age=0.000000],(73:0.000000[&&NHX:age=0.000000],((69:0.000000[&&NHX:age=0.000000],48:0.000000[&&NHX:age=0.000000])138:0.000000[&&NHX:age=0.000000],5:0.000000[&&NHX:age=0.000000])191:0.000000[&&NHX:age=0.000000])146:0.000000[&&NHX:age=0.000000])153:0.000000[&&NHX:age=0.000000])151:0.000000[&&NHX:age=0.000000])125:0.000000[&&NHX:age=0.000000])185:0.000000[&&NHX:age=0.000000])117:0.000000[&&NHX:age=0.000000])175:0.000000[&&NHX:age=0.000000])177:0.000000[&&NHX:age=0.000000],(62:0.000000[&&NHX:age=0.000000],(76:0.000000[&&NHX:age=0.000000],(56:0.000000[&&NHX:age=0.000000],((46:0.000000[&&NHX:age=0.000000],(14:0.000000[&&NHX:age=0.000000],(9:0.000000[&&NHX:age=0.000000],(77:0.000000[&&NHX:age=0.000000],(15:0.000000[&&NHX:age=0.000000],(53:0.000000[&&NHX:age=0.000000],(((87:0.000000[&&NHX:age=0.000000],(95:0.000000[&&NHX:age=0.000000],(35:0.000000[&&NHX:age=0.000000],(10:0.000000[&&NHX:age=0.000000],91:0.000000[&&NHX:age=0.000000])186:0.000000[&&NHX:age=0.000000])100:0.000000[&&NHX:age=0.000000])190:0.000000[&&NHX:age=0.000000])174:0.000000[&&NHX:age=0.000000],(57:0.000000[&&NHX:age=0.000000],29:0.000000[&&NHX:age=0.000000])114:0.000000[&&NHX:age=0.000000])144:0.000000[&&NHX:age=0.000000],(26:0.000000[&&NHX:age=0.000000],(18:0.000000[&&NHX:age=0.000000],((65:0.000000[&&NHX:age=0.000000],(42:0.000000[&&NHX:age=0.000000],(98:0.000000[&&NHX:age=0.000000],61:0.000000[&&NHX:age=0.000000])122:0.000000[&&NHX:age=0.000000])167:0.000000[&&NHX:age=0.000000])130:0.000000[&&NHX:age=0.000000],3:0.000000[&&NHX:age=0.000000])135:0.000000[&&NHX:age=0.000000])141:0.000000[&&NHX:age=0.000000])103:0.000000[&&NHX:age=0.000000])115:0.000000[&&NHX:age=0.000000])123:0.000000[&&NHX:age=0.000000])139:0.000000[&&NHX:age=0.000000])154:0.000000[&&NHX:age=0.000000])137:0.000000[&&NHX:age=0.000000])109:0.000000[&&NHX:age=0.000000])183:0.000000[&&NHX:age=0.000000],0:0.000000[&&NHX:age=0.000000])131:0.000000[&&NHX:age=0.000000])112:0.000000[&&NHX:age=0.000000])152:0.000000[&&NHX:age=0.000000])124:0.000000[&&NHX:age=0.000000])145:0.000000[&&NHX:age=0.000000])161:395.449492[&&NHX:age=0.000000])193:1959.252495[&&NHX:age=395.449492])128:0.000000[&&NHX:age=2354.701987])196:9707.106528[&&NHX:age=2354.701987])176:5982.816948[&&NHX:age=12061.808515])189:8925.972860[&&NHX:age=18044.625462])148[&&NHX:age=26970.598323];";

	//		input ="(((27:1002.805899[&&NHX:age=0.000000],(86:0.000000[&&NHX:age=0.000000],13:0.000000[&&NHX:age=0.000000])172:1002.805899[&&NHX:age=0.000000])107:0.000000[&&NHX:age=1002.805899],1:1002.805899[&&NHX:age=0.000000])101:25967.792423[&&NHX:age=1002.805899],(((75:1545.314509[&&NHX:age=0.000000],(96:1002.805899[&&NHX:age=0.000000],24:1002.805899[&&NHX:age=0.000000])192:542.508609[&&NHX:age=1002.805899])184:16499.310954[&&NHX:age=1545.314509],(22:12061.808515[&&NHX:age=0.000000],(92:12061.808515[&&NHX:age=0.000000],(74:12061.808515[&&NHX:age=0.000000],(52:8051.702368[&&NHX:age=0.000000],((88:395.449492[&&NHX:age=0.000000],4:395.449492[&&NHX:age=0.000000])193:1959.252495[&&NHX:age=395.449492],(64:122.586947[&&NHX:age=0.000000],((89:0.000000[&&NHX:age=0.000000],84:0.000000[&&NHX:age=0.000000])178:0.000000[&&NHX:age=0.000000],2:0.000000[&&NHX:age=0.000000])168:122.586947[&&NHX:age=0.000000])113:2232.115040[&&NHX:age=122.586947])128:5697.000381[&&NHX:age=2354.701987])119:4010.106147[&&NHX:age=8051.702368])176:0.000000[&&NHX:age=12061.808515])129:0.000000[&&NHX:age=12061.808515])173:5982.816948[&&NHX:age=12061.808515])189:8925.972860[&&NHX:age=18044.625462],((33:3562.255340[&&NHX:age=0.000000],(((81:0.000000[&&NHX:age=0.000000],97:0.000000[&&NHX:age=0.000000])194:0.000000[&&NHX:age=0.000000],79:0.000000[&&NHX:age=0.000000])162:395.449492[&&NHX:age=0.000000],(19:49.193481[&&NHX:age=0.000000],((17:0.000000[&&NHX:age=0.000000],47:0.000000[&&NHX:age=0.000000])104:0.000000[&&NHX:age=0.000000],12:0.000000[&&NHX:age=0.000000])187:49.193481[&&NHX:age=0.000000])185:346.256011[&&NHX:age=49.193481])149:3166.805848[&&NHX:age=395.449492])105:1801.590881[&&NHX:age=3562.255340],((54:3562.255340[&&NHX:age=0.000000],(34:49.193481[&&NHX:age=0.000000],(39:0.000000[&&NHX:age=0.000000],8:0.000000[&&NHX:age=0.000000])155:49.193481[&&NHX:age=0.000000])121:3513.061859[&&NHX:age=49.193481])108:1801.590881[&&NHX:age=3562.255340],((((94:0.000000[&&NHX:age=0.000000],((90:0.000000[&&NHX:age=0.000000],67:0.000000[&&NHX:age=0.000000])180:0.000000[&&NHX:age=0.000000],(51:0.000000[&&NHX:age=0.000000],40:0.000000[&&NHX:age=0.000000])102:0.000000[&&NHX:age=0.000000])134:0.000000[&&NHX:age=0.000000])188:0.000000[&&NHX:age=0.000000],(43:0.000000[&&NHX:age=0.000000],((((99:0.000000[&&NHX:age=0.000000],70:0.000000[&&NHX:age=0.000000])198:0.000000[&&NHX:age=0.000000],37:0.000000[&&NHX:age=0.000000])140:0.000000[&&NHX:age=0.000000],(85:0.000000[&&NHX:age=0.000000],55:0.000000[&&NHX:age=0.000000])110:0.000000[&&NHX:age=0.000000])170:0.000000[&&NHX:age=0.000000],36:0.000000[&&NHX:age=0.000000])147:0.000000[&&NHX:age=0.000000])171:0.000000[&&NHX:age=0.000000])159:0.000000[&&NHX:age=0.000000],((60:0.000000[&&NHX:age=0.000000],(63:0.000000[&&NHX:age=0.000000],(78:0.000000[&&NHX:age=0.000000],(30:0.000000[&&NHX:age=0.000000],21:0.000000[&&NHX:age=0.000000])165:0.000000[&&NHX:age=0.000000])156:0.000000[&&NHX:age=0.000000])126:0.000000[&&NHX:age=0.000000])120:0.000000[&&NHX:age=0.000000],7:0.000000[&&NHX:age=0.000000])169:0.000000[&&NHX:age=0.000000])143:1545.314509[&&NHX:age=0.000000],(11:639.178343[&&NHX:age=0.000000],(((20:0.000000[&&NHX:age=0.000000],((41:0.000000[&&NHX:age=0.000000],(71:0.000000[&&NHX:age=0.000000],(58:0.000000[&&NHX:age=0.000000],((80:0.000000[&&NHX:age=0.000000],(28:0.000000[&&NHX:age=0.000000],45:0.000000[&&NHX:age=0.000000])111:0.000000[&&NHX:age=0.000000])150:0.000000[&&NHX:age=0.000000],23:0.000000[&&NHX:age=0.000000])179:0.000000[&&NHX:age=0.000000])116:0.000000[&&NHX:age=0.000000])142:0.000000[&&NHX:age=0.000000])163:0.000000[&&NHX:age=0.000000],6:0.000000[&&NHX:age=0.000000])181:0.000000[&&NHX:age=0.000000])157:0.000000[&&NHX:age=0.000000],((49:0.000000[&&NHX:age=0.000000],44:0.000000[&&NHX:age=0.000000])195:0.000000[&&NHX:age=0.000000],((68:0.000000[&&NHX:age=0.000000],(66:0.000000[&&NHX:age=0.000000],50:0.000000[&&NHX:age=0.000000])132:0.000000[&&NHX:age=0.000000])136:0.000000[&&NHX:age=0.000000],((31:0.000000[&&NHX:age=0.000000],(25:0.000000[&&NHX:age=0.000000],16:0.000000[&&NHX:age=0.000000])197:0.000000[&&NHX:age=0.000000])106:0.000000[&&NHX:age=0.000000],(38:0.000000[&&NHX:age=0.000000],((82:0.000000[&&NHX:age=0.000000],((72:0.000000[&&NHX:age=0.000000],83:0.000000[&&NHX:age=0.000000])182:0.000000[&&NHX:age=0.000000],(59:0.000000[&&NHX:age=0.000000],93:0.000000[&&NHX:age=0.000000])118:0.000000[&&NHX:age=0.000000])166:0.000000[&&NHX:age=0.000000])164:0.000000[&&NHX:age=0.000000],(73:0.000000[&&NHX:age=0.000000],((69:0.000000[&&NHX:age=0.000000],48:0.000000[&&NHX:age=0.000000])138:0.000000[&&NHX:age=0.000000],5:0.000000[&&NHX:age=0.000000])191:0.000000[&&NHX:age=0.000000])146:0.000000[&&NHX:age=0.000000])153:0.000000[&&NHX:age=0.000000])151:0.000000[&&NHX:age=0.000000])125:0.000000[&&NHX:age=0.000000])117:0.000000[&&NHX:age=0.000000])175:0.000000[&&NHX:age=0.000000])177:0.000000[&&NHX:age=0.000000],(62:0.000000[&&NHX:age=0.000000],(76:0.000000[&&NHX:age=0.000000],(56:0.000000[&&NHX:age=0.000000],((46:0.000000[&&NHX:age=0.000000],(14:0.000000[&&NHX:age=0.000000],(9:0.000000[&&NHX:age=0.000000],(77:0.000000[&&NHX:age=0.000000],(15:0.000000[&&NHX:age=0.000000],(53:0.000000[&&NHX:age=0.000000],(((87:0.000000[&&NHX:age=0.000000],(95:0.000000[&&NHX:age=0.000000],(35:0.000000[&&NHX:age=0.000000],(10:0.000000[&&NHX:age=0.000000],91:0.000000[&&NHX:age=0.000000])186:0.000000[&&NHX:age=0.000000])100:0.000000[&&NHX:age=0.000000])190:0.000000[&&NHX:age=0.000000])174:0.000000[&&NHX:age=0.000000],(57:0.000000[&&NHX:age=0.000000],29:0.000000[&&NHX:age=0.000000])114:0.000000[&&NHX:age=0.000000])144:0.000000[&&NHX:age=0.000000],(26:0.000000[&&NHX:age=0.000000],((32:0.000000[&&NHX:age=0.000000],18:0.000000[&&NHX:age=0.000000])127:0.000000[&&NHX:age=0.000000],((65:0.000000[&&NHX:age=0.000000],(42:0.000000[&&NHX:age=0.000000],(98:0.000000[&&NHX:age=0.000000],61:0.000000[&&NHX:age=0.000000])122:0.000000[&&NHX:age=0.000000])167:0.000000[&&NHX:age=0.000000])130:0.000000[&&NHX:age=0.000000],3:0.000000[&&NHX:age=0.000000])135:0.000000[&&NHX:age=0.000000])141:0.000000[&&NHX:age=0.000000])103:0.000000[&&NHX:age=0.000000])115:0.000000[&&NHX:age=0.000000])123:0.000000[&&NHX:age=0.000000])139:0.000000[&&NHX:age=0.000000])154:0.000000[&&NHX:age=0.000000])137:0.000000[&&NHX:age=0.000000])109:0.000000[&&NHX:age=0.000000])183:0.000000[&&NHX:age=0.000000],0:0.000000[&&NHX:age=0.000000])131:0.000000[&&NHX:age=0.000000])112:0.000000[&&NHX:age=0.000000])152:0.000000[&&NHX:age=0.000000])124:0.000000[&&NHX:age=0.000000])145:639.178343[&&NHX:age=0.000000])161:906.136166[&&NHX:age=639.178343])133:3818.531712[&&NHX:age=1545.314509])196:0.000000[&&NHX:age=5363.846221])160:21606.752102[&&NHX:age=5363.846221])148:0.000000[&&NHX:age=26970.598323])158[&&NHX:age=26970.598323];";

	//		input ="(((4:5363.846221[&&NHX:age=0.000000],7:5363.846221[&&NHX:age=0.000000])13:2687.856147[&&NHX:age=5363.846221],(6:3562.255340[&&NHX:age=0.000000],1:3562.255340[&&NHX:age=0.000000])14:4489.447028[&&NHX:age=3562.255340])12:0.000000[&&NHX:age=8051.702368],(3:5363.846221[&&NHX:age=0.000000],(0:2354.701987[&&NHX:age=0.000000],(2:2354.701987[&&NHX:age=0.000000],5:2354.701987[&&NHX:age=0.000000])9:0.000000[&&NHX:age=2354.701987])10:3009.144234[&&NHX:age=2354.701987])11:2687.856147[&&NHX:age=5363.846221])8[&&NHX:age=8051.702368];";
	String[] inputArray = input.split("\\(|,|\\)");

	ArrayList<String> temp = new ArrayList<String>();

	for(int i = 0; i < inputArray.length; i++){
		if(inputArray[i].length() > 0){
			temp.add(inputArray[i]);
		}
	}

	//		System.out.println("" + temp.size());

	Double[][] age = new Double[temp.size() - 1][4];

	for(int i = 0; i < temp.size() - 1; i++){
		String[] currentThing = temp.get(i).split(":|\\[");
		String[] ageArray = currentThing[3].split("=|\\]");
		age[i][1] = Double.parseDouble(currentThing[0]); // id
		age[i][0] = Double.parseDouble(ageArray[1]); // age
		age[i][2] = Double.parseDouble(currentThing[1]); // branch length
	}

	Arrays.sort(age, new Comparator<Double[]>() {
		@Override
		public int compare(final Double[] entry1, final Double[] entry2){

			int value = entry1[0].compareTo(entry2[0]);

			if (value == 0){
				value = entry1[2].compareTo(entry2[2]);
			}

			if(value == 0){
				value = entry1[1].compareTo(entry2[1]);
			}

			return value;
		}

	});

	String[] mcraArray = temp.get(temp.size()-1).split("\\[");
	int rootBranch = Integer.parseInt(mcraArray[0]);

	int index = input.indexOf(")" + rootBranch + "[");

	Node root = fillNode(new Node(), rootBranch, index, input, age);

	Integer numTaxa = (temp.size() + 1) / 2;

	int[][] tree = new int[(temp.size() - numTaxa.intValue())][3];

	for(int i = tree.length - 1; i >= 0; i--){
		//			System.out.println(i);
		ArrayList<Node> nodeList = this.nextCoalese(root, numTaxa);
		Collections.shuffle(nodeList);

		double maxAge = -2.0;
		int ageIndex = -1;

		for(int j = 0; j < nodeList.size(); j++){
			if (nodeList.get(j).age > maxAge){
				maxAge = nodeList.get(j).age;
				ageIndex = j;
			}
		}

		Node currentNode = nodeList.get(ageIndex);
		tree[i][2] = currentNode.key;
		tree[i][0] = currentNode.lChild;
		tree[i][1] = currentNode.rChild;
		currentNode.coalesed = true;
	}

	//		for(int i = 0; i < tree.length; i++){
	//			System.out.println("" + tree[i][0] + " " + tree[i][1] + " " + tree[i][2]);
	//		}

	//		System.out.println("---");

	for(int i = 0; i < tree.length; i++){
		int currentBranch = tree[i][2];
		tree[i][2] = i + numTaxa;
		for(int j = i + 1; j < tree.length; j++){
			if(tree[j][0] == currentBranch){
				tree[j][0] = i + numTaxa;
				break;
			}
			if(tree[j][1] == currentBranch){
				tree[j][1] = i + numTaxa;
				break;
			}
		}
	}

	return tree;
}

public Node fillNode(Node currentNode, int currentBranch, int currentIndex, String input,Double age[][]){
	currentNode.key = currentBranch;

	int[] kids = findChildren(input, currentIndex);

	currentNode.lChild = kids[0];
	currentNode.rChild = kids[1];

	for(int i = 0; i < age.length; i++){
		if (age[i][1].intValue() == currentBranch){
			//				System.out.println("Found age for " + currentBranch + " " + age[i][0]);
			currentNode.age = age[i][0];
		}
	}

	int child1Index = findIndex(input, kids[0]);
	int child2Index = findIndex(input, kids[1]);

	if(kids[0] > -1){
		currentNode.leftChild = fillNode(new Node(), kids[0], child1Index, input, age);
	}

	if(kids[1] > -1){
		currentNode.rightChild = fillNode(new Node(), kids[1], child2Index, input, age);
	}

	return currentNode;
}

public int[] findChildren(String input, int index){

	int parens = 0;
	int found = 0;

	int[] children = new int[2];
	children[0] = -1;
	children[1] = -1;

	for (int j = index; j >= 0; j--){
		if (input.charAt(j) == ')'){
			parens++;
		}
		if (input.charAt(j) == '('){
			parens--;
		}
		if(parens == 1 && input.charAt(j) == ':' && input.charAt(j + 1) != 'a'){

			//				System.out.println(input.substring(j - 2, j));
			if(j > 3 && input.substring(j - 4, j).matches("\\d{4}")){
				//							System.out.println(input.substring(j - 2, j));
				if(found == 0){
					children[0] = Integer.parseInt(input.substring(j - 4, j));
				} else if(found == 1){
					children[1] = Integer.parseInt(input.substring(j - 4, j));
				}
				found++;

			} else if(j > 2 && input.substring(j - 3, j).matches("\\d{3}")){
				//							System.out.println(input.substring(j - 2, j));
				if(found == 0){
					children[0] = Integer.parseInt(input.substring(j - 3, j));
				} else if(found == 1){
					children[1] = Integer.parseInt(input.substring(j - 3, j));
				}
				found++;

			} else if(j > 1 && input.substring(j - 2, j).matches("\\d{2}")){
				//							System.out.println(input.substring(j - 2, j));
				if(found == 0){
					children[0] = Integer.parseInt(input.substring(j - 2, j));
				} else if(found == 1){
					children[1] = Integer.parseInt(input.substring(j - 2, j));
				}
				found++;

			} else if(j > 0 && input.substring(j - 1, j).matches("\\d{1}")){
				if(found == 0){
					children[0] = Integer.parseInt(input.substring(j - 1, j));
				} else if(found == 1){
					children[1] = Integer.parseInt(input.substring(j - 1, j));
				}
				found++;
			}


		}

		if(found > 1){
			break;
		}
	}

	return children;

}

public int findIndex(String input, int currentBranch){
	int result = input.indexOf(")" + currentBranch + ":");


	return result;
}

public ArrayList<Node> nextCoalese(Node currentNode, int numTaxa){
	if(currentNode.key < numTaxa){
		ArrayList<Node> returnList = new ArrayList<Node>();
		return returnList;
	}

	if(currentNode.coalesed == false){
		ArrayList<Node> returnList = new ArrayList<Node>();
		returnList.add(currentNode);
		return returnList;
	}

	ArrayList<Node> leftChild = new ArrayList<Node>();
	if(currentNode.leftChild != null){
		ArrayList<Node> temp = nextCoalese(currentNode.leftChild, numTaxa);
		leftChild.addAll(temp);
	}
	ArrayList<Node> rightChild = new ArrayList<Node>();
	if(currentNode.rightChild != null){
		ArrayList<Node> temp = nextCoalese(currentNode.rightChild, numTaxa);
		rightChild.addAll(temp);
	}

	leftChild.addAll(rightChild);

	return leftChild;

}

public double getEndPosition() {
	return endPosition;
}

public long[] getPositions() {
	return positions;
}

public ArrayList<Long> getSnpPositions() {
	return snpPositions;
}

//	public ArrayList<Integer> getSnpBranchI() {
//		return snpBranchI;
//	}

//	public ArrayList<Integer> getSnpIndexI() {
//		return snpIndexI;
//	}

public ArrayList<Double> getSnpAncestralCorrelation() {
	return snpAncestralCorrelation;
}

//	public ArrayList<Double> getSnpAncestralCorrelationI() {
//		return snpAncestralCorrelationI;
//	}

public ArrayList<Double> getSnpAncestralLambda() {
	return snpAncestralLambda;
}

//	public ArrayList<Double> getSnpAncestralLambdaI() {
//		return snpAncestralLambdaI;
//	}

public ArrayList<Double> getSnpAncestralLlr() {
	return snpAncestralLlr;
}

//	public ArrayList<Double> getSnpAncestralLlrI() {
//		return snpAncestralLlrI;
//	}

public ArrayList<Double> getSnpFreqs() {
	return snpFreqs;
}

public Double[] getBranchFreqs() {
	return branchFreqs;
}

public ArrayList<Double> getTreeMaxFreqs() {
	return treeMaxFreqs;
}

public ArrayList<Double> getSnpFreqsAtBranch() {
	return snpFreqsAtBranch;
}

public ArrayList<Integer> getSnpAncestorBranch() {
	return snpAncestorBranch;
}
}