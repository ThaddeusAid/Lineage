package uk.ac.ox.stats.aid;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.Semaphore;

public class Lineage {

	public Semaphore semaphore;

	// options
	private File optionsFile;
	private File snpFile;
	private File treeFile;
	private File legendFile;
	private File ancestorFile;
	private long[] treePositions;
	private double r2Threshold;
	private long range;
	private boolean wholeChromosome;
	private int threads;
	private File smcFile;
	private boolean usingSMC;
	private boolean usingTreefile;

	private int samples;

	// execution variables
	private int numberOfTaxa;
	private int numberOfSNPs;
	private int numberOfTrees;

	private double[][] snps;

	private long[] legendPositions;
	private String[] legendNames;
	private String[] ancestorAllele;
	private String[] alternateAllele;
	private String[] referenceSnp;

	private ArrayList<Double> snpLlrArrayList;
	private ArrayList<Double> snpLambdaArrayList;
	private ArrayList<Double> snpR2LambdaArrayList;
	private ArrayList<Integer> snpBranchesArrayList;
	private ArrayList<Double> snpCorrelationArrayList;
	private ArrayList<Integer> snpIndexArrayList;
	private ArrayList<String> snpAncestorArrayList;
	private ArrayList<String> snpNamesArrayList;
	private ArrayList<Double> snpFreqs;
	private ArrayList<Double> snpFreqsAtBranch;
	private ArrayList<Double[]> fullTreeLlrResultsArrayList;
	private ArrayList<Double[]> fullTreeLambdaResultsArrayList;
	private ArrayList<Double[]> fullTreeFreqs;
	private ArrayList<Double> treeLambda;
	private ArrayList<Integer> treeBranch;
	private ArrayList<Double> treeLlr;
	private ArrayList<Double> treeFreqs;
	private ArrayList<Long> positions;
	private TreeAnalysis[] treeAnalysisArray;

	private ArrayList<Double> snpAncestralLlrArrayList;
	private ArrayList<Double> snpAncestralLambdaArrayList;
	private ArrayList<Double> snpAncestralCorrelationArrayList;
	private ArrayList<Integer> snpAncestralBranch;

	//	private ArrayList<Double> snpLlrArrayListI;
	//	private ArrayList<Double> snpLambdaArrayListI;
	//	private ArrayList<Double> snpR2LambdaArrayListI;
	//	private ArrayList<Integer> snpBranchesArrayListI;
	//	private ArrayList<Double> snpCorrelationArrayListI;
	//	private ArrayList<Integer> snpIndexArrayListI;
	//	private ArrayList<String> snpAncestorArrayListI;
	//	private ArrayList<String> snpNamesArrayListI;
	//	private ArrayList<Double[]> fullTreeLlrResultsArrayListI;
	//	private ArrayList<Double[]> fullTreeLambdaResultsArrayListI;
	//	private ArrayList<Double> treeLambdaI;
	//	private ArrayList<Integer> treeBranchI;
	//	private ArrayList<Double> treeLlrI;
	//	private ArrayList<Double> positionsI;
	//	private TreeAnalysis[] treeAnalysisArrayI;

	private ArrayList<Double> snpAncestralLlrIArrayList;
	private ArrayList<Double> snpAncestralLambdaIArrayList;
	private ArrayList<Double> snpAncestralCorrelationIArrayList;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Lineage cmm = new Lineage(args);
		cmm.execute();
	}

	public Lineage (String[] args){
		this.optionsFile = new File("options.txt");
		//		if (args.length > 0){
		//			this.optionsFile = new File(args[0]);
		//		}
	}

	public void execute(){
		// get the settings information from the options file
		// options.txt
		// treeFile=pathtoFile
		// snpFile=pathtoFile
		// legendFile=pathtoFile
		// range=x
		// ancestorFile=pathtoFile

		long startTime = System.currentTimeMillis();

		this.snpLlrArrayList = new ArrayList<Double>();

		System.out.println("Reading Options File");

		this.treePositions = new long[0];

		this.positions = new ArrayList<Long>();
		try {
			File treeDir = new File("trees/");
			treeDir.mkdir();

			// load the options file
			this.loadOptions();

			// sort the tree positions so that it comes in order
			Arrays.sort(treePositions);

			// read the snp information into memory
			if(this.snpFile != null){
				this.loadSNPs();
			} else {
				this.snps = null;
			}

			// now that the SNPs are in memory in an array
			// load the legend file so that we have position information
			if(this.legendFile != null){
				this.loadLegend();
			}

			// load the ancestor allele information
			if(this.ancestorFile != null){
				this.loadAncestorAlleles();
			}

			System.out.println("Number of Threads: " + this.semaphore.availablePermits());

			if(usingTreefile == true){
				System.out.println("Using GENECLUSTER Tree File.");
				this.useTreefile();
			}

			if(usingSMC == true){
				System.out.println("Using ARGWeaver SMC file.");
				this.useSMCfile();
			}



		} catch (Exception e) {
			e.printStackTrace();
		}

		System.out.println("Begin Data Gathering.");

		try{
			int test = -1;
			try{
				for (int i = 0; i < this.treeAnalysisArray.length; i++){
					test = i;
					this.treeAnalysisArray[i].join();
				}
			} catch (Exception e){
				System.out.println("Nothing at " + test);
				e.printStackTrace();
			}

			this.snpLlrArrayList = new ArrayList<Double>();
			this.snpBranchesArrayList = new ArrayList<Integer>();
			this.snpLambdaArrayList = new ArrayList<Double>();
			this.snpR2LambdaArrayList = new ArrayList<Double>();
			this.snpCorrelationArrayList = new ArrayList<Double>();
			this.snpIndexArrayList = new ArrayList<Integer>();
			this.snpAncestorArrayList = new ArrayList<String>();
			this.snpNamesArrayList = new ArrayList<String>();
			this.fullTreeLlrResultsArrayList = new ArrayList<Double[]>();
			this.fullTreeLambdaResultsArrayList = new ArrayList<Double[]>();
			this.treeBranch = new ArrayList<Integer>();
			this.treeLambda = new ArrayList<Double>();
			this.treeLlr = new ArrayList<Double>();
			this.positions = new ArrayList<Long>();

			this.snpAncestralLlrArrayList = new ArrayList<Double>();
			this.snpAncestralLambdaArrayList = new ArrayList<Double>();
			this.snpAncestralCorrelationArrayList = new ArrayList<Double>();
			this.snpAncestralBranch = new ArrayList<Integer>();

			//			this.snpLlrArrayListI = new ArrayList<Double>();
			//			this.snpBranchesArrayListI = new ArrayList<Integer>();
			//			this.snpLambdaArrayListI = new ArrayList<Double>();
			//			this.snpR2LambdaArrayListI = new ArrayList<Double>();
			//			this.snpCorrelationArrayListI = new ArrayList<Double>();
			//			this.snpIndexArrayListI = new ArrayList<Integer>();
			//			this.snpAncestorArrayListI = new ArrayList<String>();
			//			this.snpNamesArrayListI = new ArrayList<String>();
			//			this.fullTreeLlrResultsArrayListI = new ArrayList<Double[]>();
			//			this.fullTreeLambdaResultsArrayListI = new ArrayList<Double[]>();
			//			this.treeBranchI = new ArrayList<Integer>();
			//			this.treeLambdaI = new ArrayList<Double>();
			//			this.treeLlrI = new ArrayList<Double>();

			this.snpAncestralLlrIArrayList = new ArrayList<Double>();
			this.snpAncestralLambdaIArrayList = new ArrayList<Double>();
			this.snpAncestralCorrelationIArrayList = new ArrayList<Double>();
			
			this.snpFreqs = new ArrayList<Double>();
			this.snpFreqsAtBranch = new ArrayList<Double>();
			this.fullTreeFreqs = new ArrayList<Double[]>();
			this.treeFreqs = new ArrayList<Double>();

			for (int i = 0; i < this.numberOfTrees; i++){
				if (this.treeAnalysisArray[i] == null){
					this.snpLlrArrayList.add(0.0);
					this.snpBranchesArrayList.add(-3);
					this.snpLambdaArrayList.add(1.0);
					this.snpR2LambdaArrayList.add(0.0);
					this.snpCorrelationArrayList.add(0.0);
					this.snpIndexArrayList.add(-1);
					this.snpAncestorArrayList.add("-");
					this.snpNamesArrayList.add("NULL");
					this.fullTreeLlrResultsArrayList.add(new Double[0]);
					this.fullTreeLambdaResultsArrayList.add(new Double[0]);
					this.treeBranch.add(-3);
					this.treeLambda.add(1.0);
					this.treeLlr.add(0.0);
					this.positions.add(0L);

					this.snpAncestralLlrArrayList.add(0.0);
					this.snpAncestralLambdaArrayList.add(1.0);
					this.snpAncestralCorrelationArrayList.add(0.0);
					this.snpAncestralBranch.add(-1);

					//					this.snpLlrArrayListI.add(0.0);
					//					this.snpBranchesArrayListI.add(-3);
					//					this.snpLambdaArrayListI.add(1.0);
					//					this.snpR2LambdaArrayListI.add(0.0);
					//					this.snpCorrelationArrayListI.add(0.0);
					//					this.snpIndexArrayListI.add(-1);
					//					this.snpAncestorArrayListI.add("-");
					//					this.snpNamesArrayListI.add("NULL");
					//					this.fullTreeLlrResultsArrayListI.add(new Double[0]);
					//					this.fullTreeLambdaResultsArrayListI.add(new Double[0]);
					//					this.treeBranchI.add(-3);
					//					this.treeLambdaI.add(1.0);
					//					this.treeLlrI.add(0.0);

					this.snpAncestralLlrIArrayList.add(0.0);
					this.snpAncestralLambdaIArrayList.add(1.0);
					this.snpAncestralCorrelationIArrayList.add(0.0);
					this.snpFreqs.add(0.0);
					this.snpFreqsAtBranch.add(0.0);
					this.treeFreqs.add(0.0);
					this.fullTreeFreqs.add(new Double[0]);
				} else {
					this.snpLlrArrayList.addAll(this.treeAnalysisArray[i].getSnpLlr());
					this.snpBranchesArrayList.addAll(this.treeAnalysisArray[i].getSnpBranch());
					this.snpLambdaArrayList.addAll(this.treeAnalysisArray[i].getSnpLambda());
					this.snpR2LambdaArrayList.addAll(this.treeAnalysisArray[i].getSnpR2());
					this.snpCorrelationArrayList.addAll(this.treeAnalysisArray[i].getSnpCorrelation());
					this.snpIndexArrayList.addAll(this.treeAnalysisArray[i].getSnpIndex());
					this.snpAncestorArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestor());
					this.snpNamesArrayList.addAll(this.treeAnalysisArray[i].getSnpName());
					this.fullTreeLlrResultsArrayList.add(this.treeAnalysisArray[i].getLlrResults());
					this.fullTreeLambdaResultsArrayList.add(this.treeAnalysisArray[i].getLambdaResults());
					this.treeBranch.addAll(this.treeAnalysisArray[i].getMaxRow());
					this.treeLambda.addAll(this.treeAnalysisArray[i].getTreeLambda());
					this.treeLlr.addAll(this.treeAnalysisArray[i].getTreeLlr());
					this.positions.addAll(this.treeAnalysisArray[i].getSnpPositions());

					this.snpAncestralLlrArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralLlr());
					this.snpAncestralLambdaArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralLambda());
					this.snpAncestralCorrelationArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralCorrelation());
					this.snpAncestralBranch.addAll(this.treeAnalysisArray[i].getSnpAncestorBranch());

					//					this.snpLlrArrayListI.addAll(this.treeAnalysisArray[i].getSnpLlrI());
					//					this.snpBranchesArrayListI.addAll(this.treeAnalysisArray[i].getSnpBranchI());
					//					this.snpLambdaArrayListI.addAll(this.treeAnalysisArray[i].getSnpLambdaI());
					//					this.snpR2LambdaArrayListI.addAll(this.treeAnalysisArray[i].getSnpR2I());
					//					this.snpCorrelationArrayListI.addAll(this.treeAnalysisArray[i].getSnpCorrelationI());
					//					this.snpIndexArrayListI.addAll(this.treeAnalysisArray[i].getSnpIndex());
					//					this.snpAncestorArrayListI.addAll(this.treeAnalysisArray[i].getSnpAncestor());
					//					this.snpNamesArrayListI.addAll(this.treeAnalysisArray[i].getSnpName());
					//										this.fullTreeLlrResultsArrayListI.add(this.treeAnalysisArray[i].getLlrResultsI());
					//										this.fullTreeLambdaResultsArrayListI.add(this.treeAnalysisArray[i].getLambdaResultsI());
					//					this.treeBranchI.addAll(this.treeAnalysisArray[i].getMaxRowI());
					//					this.treeLambdaI.addAll(this.treeAnalysisArray[i].getTreeLambdaI());
					//					this.treeLlrI.addAll(this.treeAnalysisArray[i].getTreeLlrI());

					//					this.snpAncestralLlrIArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralLlrI());
					//					this.snpAncestralLambdaIArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralLambdaI());
					//					this.snpAncestralCorrelationIArrayList.addAll(this.treeAnalysisArray[i].getSnpAncestralCorrelationI());

					this.snpFreqs.addAll(this.treeAnalysisArray[i].getSnpFreqs());
					this.snpFreqsAtBranch.addAll(this.treeAnalysisArray[i].getSnpFreqsAtBranch());
					this.treeFreqs.addAll(this.treeAnalysisArray[i].getTreeMaxFreqs());
					this.fullTreeFreqs.add(this.treeAnalysisArray[i].getBranchFreqs());
					this.treeAnalysisArray[i] = null;
				}
			}

			// writing data

			//----------------------------------------
			File datOutFile = new File("positions.dat");

			//use buffering
			OutputStream ofile = new FileOutputStream(datOutFile);
			OutputStream buffer = new BufferedOutputStream(ofile);
			ObjectOutput output = new ObjectOutputStream(buffer);

			output.writeObject(this.positions);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpLlr.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpLlrArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralLlr.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpAncestralLlrArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpBranches.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpBranchesArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpLambda.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpLambdaArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralLambda.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpAncestralLambdaArrayList);

			output.close();
			buffer.close();
			//----------------------------------------
			
			//----------------------------------------
			datOutFile = new File("snpAncestralBranches.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpAncestralBranch);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpR2.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpR2LambdaArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpCorrelations.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpCorrelationArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralCorrelations.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpAncestralCorrelationArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpIndex.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpIndexArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestry.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpAncestorArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpNames.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpNamesArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpFreqs.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpFreqs);

			output.close();
			buffer.close();
			//----------------------------------------
			
			//----------------------------------------
			datOutFile = new File("snpFreqsAtBranch.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.snpFreqsAtBranch);

			output.close();
			buffer.close();
			//----------------------------------------
			
			//----------------------------------------
			datOutFile = new File("treeMaxLlr.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.treeLlr);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxLambda.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.treeLambda);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxBranch.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.treeBranch);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxFreqs.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.treeFreqs);

			output.close();
			buffer.close();
			//----------------------------------------
			
			//----------------------------------------
			datOutFile = new File("llrFullTrees.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.fullTreeLlrResultsArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("lambdaFullTrees.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.fullTreeLambdaResultsArrayList);

			output.close();
			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("freqFullTrees.dat");

			//use buffering
			ofile = new FileOutputStream(datOutFile);
			buffer = new BufferedOutputStream(ofile);
			output = new ObjectOutputStream(buffer);

			output.writeObject(this.fullTreeFreqs);

			output.close();
			buffer.close();
			//----------------------------------------
			
			// Write Inverse Files!

			//----------------------------------------
			datOutFile = new File("snpLlrI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpLlrArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralLlrI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpAncestralLlrIArrayList);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpBranchesI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpBranchesArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpLambdaI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpLambdaArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralLambdaI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpAncestralLambdaIArrayList);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpR2I.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpR2LambdaArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpCorrelationsI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpCorrelationArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpAncestralCorrelationsI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpAncestralCorrelationIArrayList);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("snpIndexI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.snpIndexArrayListI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxLlrI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.treeLlrI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxLambdaI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.treeLambdaI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("treeMaxBranchI.dat");

			//use buffering
			//			ofile = new FileOutputStream(datOutFile);
			//			buffer = new BufferedOutputStream(ofile);
			//			output = new ObjectOutputStream(buffer);
			//
			//			output.writeObject(this.treeBranchI);
			//
			//			output.close();
			//			buffer.close();
			//----------------------------------------


			//----------------------------------------
			datOutFile = new File("llrFullTreesI.dat");

			//use buffering
			//						ofile = new FileOutputStream(datOutFile);
			//						buffer = new BufferedOutputStream(ofile);
			//						output = new ObjectOutputStream(buffer);
			//			
			//						output.writeObject(this.fullTreeLlrResultsArrayListI);
			//			
			//						output.close();
			//						buffer.close();
			//----------------------------------------

			//----------------------------------------
			datOutFile = new File("lambdaFullTreesI.dat");

			//use buffering
			//						ofile = new FileOutputStream(datOutFile);
			//						buffer = new BufferedOutputStream(ofile);
			//						output = new ObjectOutputStream(buffer);
			//			
			//						output.writeObject(this.fullTreeLambdaResultsArrayListI);
			//			
			//						output.close();
			//						buffer.close();
			//----------------------------------------

		} catch (Exception e){
			e.printStackTrace();
		}
		//		this.printAnalysis();

		System.out.println("Number of trees: " + this.numberOfTrees + " Number of mappings: " + this.snpLlrArrayList.size());
		long time = System.currentTimeMillis() - startTime;
		time /= 1000;
		long hours = time / 3600;
		long minutes = (time % 3600) / 60;
		long seconds = time % 60;

		System.out.println("Time taken: " + hours + ":" + minutes + ":" + seconds);
	}

	public void printAnalysis() {
		boolean AllSNPSPressent = true;

		for (int i = 0; i < this.numberOfTrees; i++){
			//			if (this.treeAnalysisArray[i].getBranch() == -1){
			//				AllSNPSPressent = false;
			//			}
			ArrayList<Integer> branches = this.treeAnalysisArray[i].getBranch();
			for(int j = 0; j < branches.size(); j++){
				if(branches.get(j) == -1){
					AllSNPSPressent = false;
				}
			}
		}

		File file = new File("analysis.txt");
		try {
			FileWriter fR = new FileWriter(file);
			BufferedWriter bR = new BufferedWriter(fR);

			bR.write("Analysis for " + this.treeFile.getName());
			bR.newLine();

			bR.write("All snps mapped to at least one tree: " + AllSNPSPressent);
			bR.newLine();

			bR.close();
			fR.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public int[][] makeTree(String[] tempStringArray) {


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

	public void loadOptions (){
		this.threads = 1;
		FileReader optionsFr;
		BufferedReader optionsBr;
		try {
			optionsFr = new FileReader(this.optionsFile);
			optionsBr = new BufferedReader(optionsFr);
			this.wholeChromosome = true;
			this.range = 0;

			while (optionsBr.ready()){
				String currentLine = optionsBr.readLine();
				String[] currentLineArray = currentLine.split("=");

				if (currentLineArray[0].equals("treeFile")){
					this.treeFile = new File(currentLineArray[1]);
					this.usingTreefile = true;
					this.samples = 1;
				}

				if (currentLineArray[0].equals("smcFile")){
					this.smcFile = new File(currentLineArray[1]);
					this.usingSMC = true;
					this.samples = 20;
				}

				if (currentLineArray[0].equals("snpFile")){
					this.snpFile = new File(currentLineArray[1]);
				}

				if (currentLineArray[0].equals("treePositions")){
					String[] treePositionsStringArray = currentLineArray[1].split(",");
					this.treePositions = new long[treePositionsStringArray.length];
					for (int i = 0; i < this.treePositions.length; i++){
						this.treePositions[i] = Integer.parseInt(treePositionsStringArray[i]);
					}
				}

				if (currentLineArray[0].equals("range")){
					this.range = Integer.parseInt(currentLineArray[1]); 
				}

				if (currentLineArray[0].equals("legendFile")){
					this.legendFile = new File(currentLineArray[1]);
				}

				if (currentLineArray[0].equals("r2")){
					this.r2Threshold = Double.parseDouble(currentLineArray[1]);
				}

				if (currentLineArray[0].equals("ancestorFile")){
					this.ancestorFile = new File(currentLineArray[1]);
				}



				if (currentLineArray[0].equals("wholeChromosome")){
					if (currentLineArray[1].equalsIgnoreCase("true")){
						this.wholeChromosome = true;
					} else {
						this.wholeChromosome = false;
					}
				}

				if (currentLineArray[0].equals("threads")){
					this.threads = Integer.parseInt(currentLineArray[1]);
					this.semaphore = new Semaphore(threads);
				}

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void loadSNPs(){
		System.out.println("Reading SNP File");
		FileReader snpFr;
		try {
			snpFr = new FileReader(this.snpFile);
			BufferedReader snpBr = new BufferedReader(snpFr);

			this.numberOfSNPs = 0;

			while (snpBr.ready()){
				snpBr.readLine();
				this.numberOfSNPs++;
			}

			snpBr.close();
			snpFr.close();

			snpFr = new FileReader(this.snpFile);
			snpBr = new BufferedReader(snpFr);

			this.snps = new double[0][0];

			String snpsString = snpBr.readLine();

			String[] snpsStringArray = snpsString.split(" ");

			this.numberOfTaxa = snpsStringArray.length;

			this.snps = new double[this.numberOfTaxa][this.numberOfSNPs];

			for (int i = 0; i < this.numberOfTaxa; i++){
				this.snps[i][0] = Double.parseDouble(snpsStringArray[i]);
			}

			for (int j = 1; j < this.numberOfSNPs; j++){
				if (snpBr.ready()){	
					snpsString = snpBr.readLine();

					snpsStringArray = snpsString.split(" ");

					for (int i = 0; i < this.numberOfTaxa; i++){
						this.snps[i][j] = Double.parseDouble(snpsStringArray[i]);
					}
				}
			}

			snpBr.close();
			snpFr.close();

			System.out.println("Number of haplotypes: " + this.snps.length);
			System.out.println("Number of SNPs: " + this.snps[0].length);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void loadLegend(){
		System.out.println("Reading Legend File");

		FileReader legendFr;
		try {
			legendFr = new FileReader(this.legendFile);
			BufferedReader legendBr = new BufferedReader(legendFr);

			this.legendNames = new String[this.numberOfSNPs];
			this.legendPositions = new long[this.numberOfSNPs];
			this.ancestorAllele = new String[this.numberOfSNPs];
			this.referenceSnp = new String[this.numberOfSNPs];
			this.alternateAllele = new String[this.numberOfSNPs];

			// strip the header
			legendBr.readLine();

			for (int i = 0; i < this.numberOfSNPs; i++){
				if (legendBr.ready()){
					String currentLine = legendBr.readLine().toLowerCase().trim();
					String[] legendStringArray = currentLine.split(" ");
					if(legendStringArray.length > 3){
						this.legendNames[i] = legendStringArray[0].trim();
						this.legendPositions[i] = Long.parseLong(legendStringArray[1].trim());
						this.referenceSnp[i] = legendStringArray[2].trim().toLowerCase();
						this.alternateAllele[i] = legendStringArray[3].trim().toLowerCase();
					}

				}
			}

			legendBr.close();
			legendFr.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void loadAncestorAlleles(){
		try {
			FileReader aaFr = new FileReader(this.ancestorFile);
			BufferedReader aaBr = new BufferedReader(aaFr);

			// stripping the header
			aaBr.readLine();


			for (int i = 0; i < this.legendPositions.length; i++){
				String currentLine = aaBr.readLine().toLowerCase().trim();
				String[] currentLineArray = currentLine.split("\t");
				
				if(this.legendPositions[i] != Long.parseLong(currentLineArray[1])){
					System.out.println("Legend != Anc");
				}
				
				if(currentLineArray.length > 3){
					this.ancestorAllele[i] = currentLineArray[4].trim();
				}
			}

			aaBr.close();
			aaFr.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void useTreefile(){
		try{
			// now we need the tree file and to start the processing
			FileReader treeFr = new FileReader(this.treeFile);
			BufferedReader treeBr = new BufferedReader(treeFr);

			while (treeBr.ready()){
				treeBr.readLine();
				this.numberOfTrees++;
			}

			System.out.println("Number of Trees: " + this.numberOfTrees);

			treeBr.close();
			treeFr.close();

			this.treeAnalysisArray = new TreeAnalysis[this.numberOfTrees];

			treeFr = new FileReader(this.treeFile);
			treeBr = new BufferedReader(treeFr);

			//Get the jvm heap size.
			long heapSize = Runtime.getRuntime().totalMemory();

			//Print the jvm heap size.
			System.out.println("Heap Size = " + heapSize);

			System.out.println("Starting to analyze trees");

			for (int treeIndex = 0; treeIndex < this.treeAnalysisArray.length; treeIndex++){
				String currentTreeString = treeBr.readLine();
				String[] currentTreeStringArray = currentTreeString.split(";");
				String[] treeHeader = currentTreeStringArray[0].split(" ");

				double i = Double.parseDouble(treeHeader[0]);

				if (this.wholeChromosome == false){
					boolean skip = true;

					for (int j = 0; j < this.treePositions.length; j++){
						if (i == this.treePositions[j]){
							skip = false;
						}
					}

					if (skip == true){
						continue;
					}
				}

				this.treeAnalysisArray[treeIndex] = new TreeAnalysis(snps, currentTreeStringArray, legendPositions, numberOfTaxa, range, (long) i, r2Threshold, this.semaphore, this.referenceSnp, this.alternateAllele, this.ancestorAllele, this.legendNames);
				try {
					this.treeAnalysisArray[treeIndex].semaphore.acquire();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				this.treeAnalysisArray[treeIndex].start();

				//System.out.println("" + i);
				this.positions.add((long) i);

			}

			
			for(int i =0; i < this.treeAnalysisArray.length; i++){
				this.treeAnalysisArray[i].join();
			}
			
			treeBr.close();
			treeFr.close();
			

			
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void useSMCfile(){
		try{
			// now we need the tree file and to start the processing
			FileReader smcFr = new FileReader(this.smcFile);
			BufferedReader smcBr = new BufferedReader(smcFr);

			while (smcBr.ready()){
				String currentLine = smcBr.readLine();
				if(currentLine.startsWith("TREE")){
					this.numberOfTrees++;
				}
			}

			System.out.println("Number of Trees: " + this.numberOfTrees);

			smcBr.close();
			smcFr.close();

			this.treeAnalysisArray = new TreeAnalysis[this.numberOfTrees];

			smcFr = new FileReader(this.smcFile);
			smcBr = new BufferedReader(smcFr);

			//Get the jvm heap size.
			long heapSize = Runtime.getRuntime().totalMemory();

			//Print the jvm heap size.
			System.out.println("Heap Size = " + heapSize);

			System.out.println("Starting to analyze trees");

			String members = smcBr.readLine();
			String[] membersArray = members.split("\t");

			int[] membersIntArray = new int[membersArray.length - 1];

			for(int i = 1; i < membersArray.length; i++){
				membersIntArray[i - 1] = Integer.parseInt(membersArray[i]);
			}

			for (int treeIndex = 0; treeIndex < this.treeAnalysisArray.length; treeIndex++){
				String currentsmcString = smcBr.readLine();

				while(!currentsmcString.startsWith("TREE")){
					currentsmcString = smcBr.readLine();
				}

				String[] currentsmcStringArray = currentsmcString.split("\t");
				String treeHeader = currentsmcStringArray[1];

				double i = Double.parseDouble(treeHeader);

				if (this.wholeChromosome == false){
					boolean skip = true;

					for (int j = 0; j < this.treePositions.length; j++){
						if (i == this.treePositions[j]){
							skip = false;
						}
					}

					if (skip == true){
						continue;
					}
				}

				this.treeAnalysisArray[treeIndex] = new TreeAnalysis(snps, currentsmcStringArray, legendPositions, numberOfTaxa, range, (long) i, r2Threshold, this.semaphore, this.referenceSnp, this.alternateAllele, this.ancestorAllele, this.legendNames, membersIntArray);
				try {
					this.treeAnalysisArray[treeIndex].semaphore.acquire();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				this.treeAnalysisArray[treeIndex].start();

				//System.out.println("" + i);
				this.positions.add((long) i);

			}

			smcBr.close();
			smcFr.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
