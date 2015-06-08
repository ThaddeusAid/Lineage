package uk.ac.ox.stats.aid;

import java.util.Arrays;

public class AidZhantikinMyers extends Thread {

	private int[][] tree;
	private int size;
	private int treeSize;
	private double lambdaLimit;
	private double lambdaStep;
	private int[][] matrix3;
	private int lambdaCols;
	private double maxScoreSelection;
	private int maxScoreRow;
	private int[][] matrix2;
	private double maxScoreNoSelection;
	private double bigLambda;
	private Double[] llrResults;
	private Double[] lambdaResults;
	private Double[] branchFreqs;

	public AidZhantikinMyers(int[][] inputTree, int inputSize, int inputTreeSize, double mL, double ls) {
        this.tree = inputTree;
        this.size = inputSize;
        this.treeSize = inputTreeSize;
        this.lambdaLimit = mL;
        this.lambdaStep = ls;
    }

	@Override
	public synchronized void start() {
		super.start();
	}
	
	@Override
	public void run() {
		super.run();
		
		this.analyse();
	}
	
	public void analyse() {

		this.lambdaCols = (int)(((this.lambdaLimit - 1.0) / this.lambdaStep) + 1.0);
        this.makeMatrix2();
        this.makeMatrix3();
        this.branchFreqs = new Double[this.matrix3.length];
        
        for(int i = 0; i < this.branchFreqs.length; i++){
        	this.branchFreqs[i] = (double)this.matrix3[i][0] / (double)this.size;
        }

        // making the hit_matrix_B like in Zalmat's code
        int[] hit_matrix_B = new int[this.size];
        
        System.arraycopy(this.matrix3[this.treeSize - 1], 0, hit_matrix_B, 0, this.size);
        
        // now that we have matrix 3 we can find out the log likelihood ratio statistic
        double[][] pArray = new double[this.treeSize][this.lambdaCols];
        
        // set the array to negative infinity
        for (int i = 0; i < this.treeSize; i++){
            Arrays.fill(pArray[i], Double.NEGATIVE_INFINITY);
        }
        
        double[] value2 = new double[this.treeSize];
        int column = 0;

        double part3 = 0.0;
        
        for (int i = 1; i < this.size; i++) {
            part3 += Math.log(i);
        }

        for (double lambda = 1.0; lambda <= this.lambdaLimit; lambda += this.lambdaStep) {

            for (int row = 0; row < this.treeSize; row++) {
                int N = this.matrix3[row][0];

                // equation is (N - 2) * log(lambda) - sum of (log (i - li + li*lambda)) plus sum of (log(i))
                double part1 = (N - 2) * Math.log(lambda);  // log(lamba^(N-2))

                double part2 = 0.0;

                for (int i = 1; i < this.treeSize; i++) {
                    part2 += Math.log((double) hit_matrix_B[i] + (double) this.matrix3[row][i] * (lambda - 1));
                }

                
                
                double value = part1 - part2 + part3;

                // add them all together

                pArray[row][column] = value;

                if (lambda == 1.0) {
                    value2[row] = value;
                }
            }

            column++;
        }

        // find max score
        this.maxScoreSelection = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < this.size - 1; i++) {
            for (int j = 0; j < this.lambdaCols; j++) {
                if (pArray[i][j] >= maxScoreSelection) {
                    maxScoreSelection = pArray[i][j];
                    this.maxScoreRow = i;
                }
            }
        }
        
        maxScoreNoSelection = value2[maxScoreRow];
//        System.out.println("hit_matrix_B[0] = " + hit_matrix_B[0] + "Value1 = " + maxScoreSelection + " Value2 = " + value2[maxScoreRow] + "BigLambda = " + (2.0 * (maxScoreSelection - maxScoreNoSelection)));
        this.bigLambda = maxScoreSelection - maxScoreNoSelection;

        this.bigLambda *= 2;
        
        // collect the max llr for each branch
        this.llrResults = new Double[this.size - 1];
        this.lambdaResults = new Double[this.size -1];
        for (int i = 0; i < this.size - 1; i++){
        	double limit = Double.NEGATIVE_INFINITY;
        	int index = Integer.MIN_VALUE;
        	for (int j = 0; j < this.lambdaCols; j++){
        		if (pArray[i][j] >= limit){
        			limit = pArray[i][j];
        			index = j;
        		}
        	}
        	limit -= pArray[i][0];
        	limit *= 2;
        	this.llrResults[i] = limit;
        	this.lambdaResults[i] = 1.0 + (this.lambdaStep * (double)index);
        }
    }

    public void makeMatrix2(){
        int matrix2rows = 2 * this.size - 1;
        
//        System.out.println("Matrix Rows: " + matrix2rows + " this.size: " + this.size);
        this.matrix2 = new int[matrix2rows][this.size];
        
        // set first column of lineages to 1
        for (int i = 0; i < this.size; i++){
            this.matrix2[i][0] = 1;
        }
        
        // grow the life of the lineages
        for (int i = 0; i < this.size; i++){
            //starting at column 1 and tree row 0 (column - 1)
            for (int j = 1; j < this.size; j++){
                if (this.tree[j - 1][0] == i || this.tree[j-1][1] == i){
                    break;
                } else {
                    this.matrix2[i][j] = 1;
                }
            }
        }
        
        // fill the bottom part of the matrix
        for (int i = 0; i < this.treeSize; i++){
            int m = this.tree[i][2];
            int n = this.tree[i][0];
            int o = this.tree[i][1];
            
            for (int j = 0; j < this.size; j++){
                this.matrix2[m][j] = this.matrix2[n][j] + this.matrix2[o][j];
            }
            
            for (int j = i; j < this.treeSize; j++){
                if (this.tree[j][0] == m || this.tree[j][1] == m){
                    break;
                } else {
                    this.matrix2[m][j + 1] = 1;
                }
            }
        }
    }
    
    public void makeMatrix3(){
        int matrix2Rows = 2 * this.size - 1;
        this.matrix3 = new int[this.treeSize][this.size];
        
        // matrix 2 is filled in and need to be transferred to matrix3
        for (int i = this.size; i < matrix2Rows; i++) {
            System.arraycopy(this.matrix2[i], 0, this.matrix3[i - this.size], 0, this.size);
        }

        // remove the extra 1s and 2s from matrix3
        for (int i = 0; i < this.size - 1; i++) {
            for (int j = 0; j < this.size; j++) {
                if (this.matrix3[i][j] < 2){
                	this.matrix3[i][j] = 0;
                }
            }
        }
    }
    
    public void makeTree(String[] tempStringArray) {

//                System.out.println("SNP " + tempStringArray[0]);
        this.size = Integer.parseInt(tempStringArray[1].substring(1));
        tempStringArray = tempStringArray[2].split(",");
        this.treeSize = tempStringArray.length;

        int[][] currentTree = new int[tempStringArray.length][3];

        for (int j = 0; j < tempStringArray.length; j++) {
            String[] branches = tempStringArray[j].split(" ");
            currentTree[j][0] = Integer.parseInt(branches[1]) - 1;
            currentTree[j][1] = Integer.parseInt(branches[2]) - 1;
            currentTree[j][2] = this.size + j;
//                System.out.println(this.treesArray[i][j][0] + "," + this.treesArray[i][j][1] + " " + this.treesArray[i][j][2]);
        }

        this.tree = currentTree;
    }
    
    public void makeTree(String tempString) {
        String[] tempStringArray = tempString.split(";");

//                System.out.println("SNP " + tempStringArray[0]);
        this.size = Integer.parseInt(tempStringArray[1].substring(1));
        tempStringArray = tempStringArray[2].split(",");
        this.treeSize = tempStringArray.length;

        int[][] currentTree = new int[tempStringArray.length][3];

        for (int j = 0; j < tempStringArray.length; j++) {
            String[] branches = tempStringArray[j].split(" ");
            currentTree[j][0] = Integer.parseInt(branches[1]) - 1;
            currentTree[j][1] = Integer.parseInt(branches[2]) - 1;
            currentTree[j][2] = this.size + j;
//                System.out.println(this.treesArray[i][j][0] + "," + this.treesArray[i][j][1] + " " + this.treesArray[i][j][2]);
        }

        this.tree = currentTree;
    }
    
    public double getBigLambda() {
		return bigLambda;
	}
    
    public Double[] getLlrResults() {
		return llrResults;
	}
    
    public Double[] getLambdaResults() {
		return lambdaResults;
	}
    
    public Double[] getBranchFreqs() {
		return branchFreqs;
	}
}
