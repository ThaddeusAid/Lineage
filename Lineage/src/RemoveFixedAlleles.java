import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

// usage java RemoveFixedAlleles pre=ceuChr6 legendFile=ceuChr6.impute.legend hapFile=ceuChr6.impute.hap

public class RemoveFixedAlleles {

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		String hap = "";
		String legend = "";
		String anc = "";
		String pre = "";
		
		for (int i = 0; i < args.length; i++){
			String[] currentArgArray = args[i].split("=");
			if (currentArgArray[0].contains("hapFile")){
				hap = currentArgArray[1].trim();
			}
			
			if (currentArgArray[0].contains("legendFile")){
				legend = currentArgArray[1].trim();
			}
			
			if (currentArgArray[0].contains("pre")){
				pre = currentArgArray[1].trim();
			}
			
			if (currentArgArray[0].contains("ancestorFile")){
				anc = currentArgArray[1].trim();
			}
		}
		
		File legendFile = new File(legend);
		File hapFile = new File(hap);
		File ancestorFile = new File(anc);
		
		try{
			FileReader legendFr = new FileReader(legendFile);
			BufferedReader legendBr = new BufferedReader(legendFr);
			
			FileReader hapFr = new FileReader(hapFile);
			BufferedReader hapBr = new BufferedReader(hapFr);
			
			FileReader ancFr = new FileReader(ancestorFile);
			BufferedReader ancBr = new BufferedReader(ancFr);
			
			FileWriter legendFw = new FileWriter(new File(pre + ".thinned.legend"));
			BufferedWriter legendBw = new BufferedWriter(legendFw);
			
			FileWriter hapFw = new FileWriter(new File(pre + ".thinned.hap"));
			BufferedWriter hapBw = new BufferedWriter(hapFw);
			
			FileWriter ancFw = new FileWriter(new File(pre + ".thinned.anc"));
			BufferedWriter ancBw = new BufferedWriter(ancFw);
			
			String temp = legendBr.readLine().trim();
			
			legendBw.write(temp);
			legendBw.newLine();
			
			temp = ancBr.readLine().trim();
			
			ancBw.write(temp);
			ancBw.newLine();
			
			while(legendBr.ready() && hapBr.ready() && ancBr.ready()){
				String hapTemp = hapBr.readLine();
				String legendTemp = legendBr.readLine();
				String ancTemp = ancBr.readLine();
				
				if (hapTemp.contains("0") && hapTemp.contains("1")){
					String[] hapArray = hapTemp.split(" ");
					
					int total = 0;
					
					for(int i = 0; i < hapArray.length; i++){
						total += Integer.parseInt(hapArray[i]);
					}
					
					String[] ancArray = ancTemp.split("\t");
					
					if (ancArray[2].equalsIgnoreCase(ancArray[4].split("|")[0])){
						if(total <= 1){
							continue;
						}
					} else if (ancArray[3].equalsIgnoreCase(ancArray[4].split("|")[0])){
						if(total >= hapArray.length - 1){
							continue;
						}
					}
					
					legendBw.write(legendTemp);
					legendBw.newLine();
					hapBw.write(hapTemp);
					hapBw.newLine();
					ancBw.write(ancTemp);
					ancBw.newLine();
				}
			}
			
			hapBw.close();
			hapFw.close();
			
			hapBr.close();
			hapFr.close();
			
			legendBw.close();
			legendFw.close();
			
			legendBr.close();
			legendFr.close();
			
			ancBw.close();
			ancFw.close();
			
			ancBr.close();
			ancFr.close();
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}

}
