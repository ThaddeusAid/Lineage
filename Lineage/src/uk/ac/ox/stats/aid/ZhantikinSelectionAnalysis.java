/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ox.stats.aid;

import java.util.Arrays;

/**
 *
 * @author Thad
 */
public class ZhantikinSelectionAnalysis extends Thread {

    private int size;
    private int treeSize;
    private int[][] tree;
    private int[][] matrix2;
    private int[][] matrix3;
    private double bigLambda;
    private int maxScoreRow;
    private int maxScoreColumn;
    private double maxScoreSelection;
    private double maxScoreNoSelection;
    private double lambdaLimit;
    private int lambdaCols;
    private double lambdaStep;
    private Double[] llrResults;
    private Double[] lambdaResults;

    public ZhantikinSelectionAnalysis() {
    }

    public ZhantikinSelectionAnalysis(int[][] inputTree, int inputSize, int inputTreeSize, double mL, double ls) {
        this.tree = inputTree;
        this.size = inputSize;
        this.treeSize = inputTreeSize;
        this.bigLambda = mL;
        this.lambdaStep = ls;
    }

    public void analyzeSingle() {

        this.makeMatrix2();
        this.makeMatrix3();

        // making the hit_matrix_B like in Zalmat's code
        int[] hit_matrix_B = new int[this.size];
        
        System.arraycopy(this.matrix3[this.treeSize - 1], 0, hit_matrix_B, 0, this.size);

        // now that we have matrix 3 we can find out the log likelihood ratio statistic
        double[][] pArray = new double[this.treeSize][this.lambdaCols];
        double[] value2 = new double[this.treeSize];
        int column = 0;

        double part3 = 0.0;

        for (int i = 0; i < this.size; i++) {
            part3 += Math.log(i + 1);
        }

        for (double lambda = 1.0; lambda <= this.lambdaLimit; lambda += this.lambdaStep) {

            for (int row = 0; row < this.treeSize; row++) {
                int N = this.matrix3[row][0];

                // equation is (N - 1) * log(lambda) - sum of (log (i - li + li*lambda)) plus sum of (log(i))
                double part1 = (N - 1) * Math.log(lambda);

                double part2 = 0.0;

                for (int i = 1; i < this.size; i++) {
                    part2 += Math.log((double) hit_matrix_B[i] - (double)this.matrix3[row][i] + (double) this.matrix3[row][i] * lambda);
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
        this.maxScoreRow = -1;

        for (int i = 0; i < this.size - 1; i++) {
            for (int j = 0; j < this.lambdaCols; j++) {
                if (pArray[i][j] > maxScoreSelection) {
                    maxScoreSelection = pArray[i][j];
                    maxScoreRow = i;
                    maxScoreColumn = j;
                }
            }
        }

        maxScoreNoSelection = value2[maxScoreRow];
//        System.out.println("hit_matrix_B[0] = " + hit_matrix_B[0] + "Value1 = " + maxScoreSelection + " Value2 = " + value2[maxScoreRow] + "BigLambda = " + (2.0 * (maxScoreSelection - maxScoreNoSelection)));
        this.bigLambda = maxScoreSelection - maxScoreNoSelection;

        this.bigLambda *= 2;
    }

    public int[][] getMatrix2() {
        return matrix2;
    }

    public int[][] getMatrix3() {
        return matrix3;
    }

    public int getSize() {
        return size;
    }

    public int[][] getTree() {
        return tree;
    }

    public int getTreeSize() {
        return treeSize;
    }

    public void setMatrix2(int[][] matrix2) {
        this.matrix2 = matrix2;
    }

    public void setMatrix3(int[][] matrix3) {
        this.matrix3 = matrix3;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public void setTree(int[][] tree) {
        this.tree = tree;
    }

    public void setTreeSize(int treeSize) {
        this.treeSize = treeSize;
    }

    public double getMaxLambda() {
        return bigLambda;
    }

    public void setMaxLambda(double maxLambda) {
        this.bigLambda = maxLambda;
    }

    public double getBigLambda() {
		return bigLambda;
	}
    
    public void setBigLambda(double bigLambda) {
		this.bigLambda = bigLambda;
	}
    
    public synchronized void start(double lambda, double lambdastep) {
        super.start();

        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
    }
    
    public synchronized void run(){
    	super.run();
    	this.analyzeSingle2();
    }
    
    @Override
    public synchronized void start() {
    	// TODO Auto-generated method stub
    	super.start();
    }
    
    public synchronized void run(double lambda, double lambdastep){
        super.run();
        
        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
        this.analyzeSingle2();
    }
    
    public synchronized void start(String[] treeStringArray, double lambda, double lambdastep){
        super.start();
        
        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
        this.makeTree(treeStringArray);
        this.analyzeSingle2();
    }
    
    public synchronized void run(String[] treeStringArray, double lambda, double lambdastep){
        super.run();
        
        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
        this.makeTree(treeStringArray);
        this.analyzeSingle2();
    }
    
    public synchronized void start(String treeString, double lambda, double lambdastep){
        super.start();
        
        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
        this.makeTree(treeString);
        this.analyzeSingle2();
    }
    
    public synchronized void run(String treeString, double lambda, double lambdastep){
        super.run();
        
        this.lambdaLimit = lambda;
        this.lambdaStep = lambdastep;
        
        this.lambdaCols = (int)((this.lambdaLimit - 1) / this.lambdaStep) + 1;
        
        this.makeTree(treeString);
        this.analyzeSingle2();
    }
    
    
    
    public void analyzeSingle2() {

        this.makeMatrix2();
        this.makeMatrix3();

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

        for (int i = 0; i < this.size; i++) {
            part3 += Math.log(i + 1);
        }

        for (double lambda = 1.0; lambda <= this.lambdaLimit; lambda += this.lambdaStep) {

            for (int row = 0; row < this.treeSize; row++) {
                int N = this.matrix3[row][0];

                // equation is (N - 2) * log(lambda) - sum of (log (i - li + li*lambda)) plus sum of (log(i))
                double part1 = (N - 2) * Math.log(lambda);

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
                    this.maxScoreColumn = j;
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

        // remove the extra 1s from matrix3
        for (int i = 0; i < this.size - 1; i++) {
            for (int j = i + 2; j < this.size; j++) {
                this.matrix3[i][j] = 0;
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

    public int getMaxScoreColumn() {
        return maxScoreColumn;
    }

    public double getMaxScoreNoSelection() {
        return maxScoreNoSelection;
    }

    public int getMaxScoreRow() {
        return maxScoreRow;
    }

    public double getMaxScoreSelection() {
        return maxScoreSelection;
    }

    public void setMaxScoreColumn(int maxScoreColumn) {
        this.maxScoreColumn = maxScoreColumn;
    }

    public void setMaxScoreNoSelection(double maxScoreNoSelection) {
        this.maxScoreNoSelection = maxScoreNoSelection;
    }

    public void setMaxScoreRow(int maxScoreRow) {
        this.maxScoreRow = maxScoreRow;
    }

    public void setMaxScoreSelection(double maxScoreSelection) {
        this.maxScoreSelection = maxScoreSelection;
    }
    
    public Double[] getLlrResults() {
		return llrResults;
	}
    
    public void setLlrResults(Double[] llrResults) {
		this.llrResults = llrResults;
	}
    
    public int getLambdaCols() {
		return lambdaCols;
	}
    public double getLambdaLimit() {
		return lambdaLimit;
	}
    public Double[] getLambdaResults() {
		return lambdaResults;
	}
    public double getLambdaStep() {
		return lambdaStep;
	}
}
