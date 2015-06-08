import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class FixAncestralInformation {
	
	private File legendFile;
	private File snpFile;
	private File ancestryFile;

	public static void main(String[] args) {
		FixAncestralInformation fai = new FixAncestralInformation();
		fai.execute();
	}
	
	public FixAncestralInformation(){
		
	}

	public void execute(){
		
		try {
			// load options
			File optionsFile = new File("options.txt");
			
			if(optionsFile.exists()){
				System.out.println("options.txt exists!");
			}
			FileReader optionsFr = new FileReader(optionsFile);
			BufferedReader optionsBr = new BufferedReader(optionsFr);
			
			while(optionsBr.ready()){
				String currentLine = optionsBr.readLine();
				String[] currentLineArray = currentLine.split("=");
				
				if (currentLineArray[0].contentEquals("legendFile")){
					this.legendFile = new File(currentLineArray[1].trim());
					System.out.println("Legend File: " + this.legendFile.getAbsolutePath());
				}
				
				if (currentLineArray[0].contentEquals("snpFile")){
					this.snpFile = new File(currentLineArray[1].trim());
					System.out.println("SNP file: " + this.snpFile.getAbsoluteFile());
				}
				
				if (currentLineArray[0].contentEquals("ancestryFile")){
					this.ancestryFile = new File(currentLineArray[1].trim());
					System.out.println("Ancestry File: " + this.ancestryFile.getAbsolutePath());
					
				}
			}
			
			optionsBr.close();
			optionsFr.close();
			
			// open legend file and snp file
			// legend file - name pos ref alt
			// snp file - snps on rows sep = " "
			// ancestry file - chr pos ref alt anc
			int linesRemoved = 0;
			int linesKept = 0;
			int linesAltered = 0;
			
			ArrayList<Double> frequency = new ArrayList<Double>();
			ArrayList<String> positions = new ArrayList<String>();
			
			FileReader ancestryFr = new FileReader(ancestryFile);
			BufferedReader ancestryBr = new BufferedReader(ancestryFr);
			
			FileReader legendFr = new FileReader(legendFile);
			BufferedReader legendBr = new BufferedReader(legendFr);
			
			FileReader snpFr = new FileReader(snpFile);
			BufferedReader snpBr = new BufferedReader(snpFr);
			
			FileWriter ancestryFw = new FileWriter(new File(this.ancestryFile.getAbsoluteFile() + ".ancestry"));
			BufferedWriter ancestryBw = new BufferedWriter(ancestryFw);
			
			FileWriter legendFw = new FileWriter(legendFile.getAbsoluteFile() + ".ancestry");
			BufferedWriter legendBw = new BufferedWriter(legendFw);
			
			FileWriter snpFw = new FileWriter(snpFile.getAbsoluteFile() + ".ancestry");
			BufferedWriter snpBw = new BufferedWriter(snpFw);
			
			//strip headers
			ancestryBr.readLine();
			legendBr.readLine();
			
			while(ancestryBr.ready() && legendBr.ready() && snpBr.ready()){
				String currentAncestryLine = ancestryBr.readLine().trim();
				String currentLegendLine = legendBr.readLine().trim();
				String currentSnpLine = snpBr.readLine().trim();
				
				String[] currentAncestryLineArray = currentAncestryLine.split("\t");
				
				if(currentAncestryLineArray[2].length() > 1 || 
						currentAncestryLineArray[3].length() > 1 ||
						currentAncestryLineArray[4].contains("-") ||
						currentAncestryLineArray[4].contains(".") ||
						currentAncestryLineArray[4].contains("n") ||
						!currentSnpLine.contains("1") ||
						!currentSnpLine.contains("0")){
					linesRemoved += 1;
					continue;
				}
				
				if(!currentAncestryLineArray[2].contentEquals(currentAncestryLineArray[4]) &&
						!currentAncestryLineArray[3].contentEquals(currentAncestryLineArray[4])){
					linesRemoved += 1;
					continue;
				}
				
				if (currentAncestryLineArray[2].contentEquals(currentAncestryLineArray[4])){

					int ones = 0;
					int zeros = 0;
					
					for (int i = 0; i < currentSnpLine.length(); i++){
						if (currentSnpLine.charAt(i) == '0'){
							zeros += 1;
						}
						
						if (currentSnpLine.charAt(i) == '1'){
							ones += 1;
						}
					}
					
//					if (ones < 2){
//						linesRemoved +=1;
//						continue;
//					}
					positions.add(currentAncestryLineArray[1]);
					frequency.add((double)ones / (double)(ones + zeros));
					
					ancestryBw.write(currentAncestryLine);
					ancestryBw.newLine();
					
					legendBw.write(currentLegendLine);
					legendBw.newLine();
					
					snpBw.write(currentSnpLine);
					snpBw.newLine();
					
					linesKept += 1;
				}
				
				if (currentAncestryLineArray[3].contentEquals(currentAncestryLineArray[4])){
					currentSnpLine.replaceAll("1", "2");
					currentSnpLine.replaceAll("0", "1");
					currentSnpLine.replaceAll("2", "0");
					
					int ones = 0;
					int zeros = 0;
					
					for (int i = 0; i < currentSnpLine.length(); i++){
						if (currentSnpLine.charAt(i) == '0'){
							zeros += 1;
						}
						
						if (currentSnpLine.charAt(i) == '1'){
							ones += 1;
						}
					}
					
//					if(ones < 2){
//						linesRemoved += 1;
//						continue;
//					}
					
					positions.add(currentAncestryLineArray[1]);
					frequency.add((double)ones / (double)(ones + zeros));
					
					ancestryBw.write(currentAncestryLineArray[0] + " " + 
							currentAncestryLineArray[1] + "\t" +
							currentAncestryLineArray[3] + "\t" +
							currentAncestryLineArray[2] + "\t" +
							currentAncestryLineArray[4]);
					ancestryBw.newLine();
					
					String[] currentLegendLineArray = currentLegendLine.split(" ");
					legendBw.write(currentLegendLineArray[0] + " " +
							currentLegendLineArray[1] + " " +
							currentLegendLineArray[3] + " " +
							currentLegendLineArray[2]);
					legendBw.newLine();
					
					snpBw.write(currentSnpLine);
					snpBw.newLine();
					
					linesAltered += 1;
				}
			}
			
			snpBw.close();
			snpFw.close();
			
			snpBr.close();
			snpFr.close();
			
			legendBw.close();
			legendFw.close();
			
			legendBr.close();
			legendFr.close();
			
			ancestryBw.close();
			ancestryFw.close();
			
			ancestryBr.close();
			ancestryFr.close();
			
			// write analysis
			FileWriter analysisFw = new FileWriter(new File("analysis.txt"));
			BufferedWriter analysisBw = new BufferedWriter(analysisFw);
			
			analysisBw.write("Lines Altered: " + linesAltered);
			analysisBw.newLine();
			
			analysisBw.write("Lines Kept: " + linesKept);
			analysisBw.newLine();
			
			analysisBw.write("Lines Removed: " + linesRemoved);
			analysisBw.newLine();
			
			analysisBw.close();
			analysisFw.close();
			
			// write frequencies
			FileWriter freqFw = new FileWriter(new File("frequencies.txt"));
			BufferedWriter freqBw = new BufferedWriter(freqFw);
			
			for(int i = 0; i < frequency.size(); i++){
				DecimalFormat df = new DecimalFormat("0.000000");
				freqBw.write(positions.get(i) + "," + df.format(frequency.get(i).doubleValue()));
				freqBw.newLine();
			}
			
			freqBw.close();
			freqFw.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		
		
		
	}
}
