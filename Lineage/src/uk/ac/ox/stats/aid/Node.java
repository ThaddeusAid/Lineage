package uk.ac.ox.stats.aid;
class Node {

	int key;
	String name;
	double age;
	int lChild;
	int rChild;
	boolean coalesed;

	Node leftChild;
	Node rightChild;

	Node(){
		key = -1;
		name = "";
		age = -0.5;
		lChild = -1;
		rChild = -1;
		leftChild = null;
		rightChild = null;
		coalesed = false;
	}

	Node(int key, Double a, int lChild, int rChild) {
		this.key = key;
		this.age = a;
		this.rChild = rChild;
		this.lChild = lChild;

	}

	Node(int key, String name){
		this.key = key;
		this.name = name;
	}

	public String toString() {	 	        
		return name + " has the key " + key;

		/*
		 * return name + " has the key " + key + "\nLeft Child: " + leftChild +         * "\nRight Child: " + rightChild + "\n";
		 */

	}

	

}