import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class CorrectByWT3 {
	public static void main(String[] args) throws IOException{
		//Proteomic
		File f1 = new File("D:\\Dmelanogaster\\Proteome\\new_data1.txt");
		File f2 = new File("D:\\Dmelanogaster\\Proteome\\new_data2.txt");
		if(f1.exists() && f2.exists()){
			BufferedReader br1 = new BufferedReader(new FileReader(f1));
			BufferedReader br2 = new BufferedReader(new FileReader(f2));
			
			BufferedWriter bw1 = new BufferedWriter(new FileWriter
					("D:\\Dmelanogaster\\Proteome\\new_data1_correct.txt"));
			String s1;
			Map<String,Double> hm2yw3 = new HashMap<String,Double>();
			while((s1 = br2.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String id = sp1[0];
				double yw3 = Double.parseDouble(sp1[1]);
				if(yw3 != 0.0){
					hm2yw3.put(id, yw3);
				}
				
			}
			
			while((s1 = br1.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String id = sp1[0];
				if(!hm2yw3.containsKey(id)){
					continue;
				}else{
					double yw3 = Double.parseDouble(sp1[2]);
					if(yw3 == 0.0){
						continue;
					}
					double foldchange = hm2yw3.get(id)/yw3;
					bw1.write(id);
					for(int i=1;i<sp1.length;i++){
					if(i == 2){
						continue;
					}
						double inten = Double.parseDouble(sp1[i]);
						double newinten = inten * foldchange;
						bw1.write("\t" + String.format("%.2f", newinten));
					}
					bw1.newLine();
				}
			}
			
			bw1.flush();
			bw1.close();
		}
		//phosphoproteomic
		/*File f1 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data1.txt");
		File f2 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data2.txt");
		if(f1.exists() && f2.exists()){
			BufferedReader br1 = new BufferedReader(new FileReader(f1));
			BufferedReader br2 = new BufferedReader(new FileReader(f2));
			
			BufferedWriter bw1 = new BufferedWriter(new FileWriter
					("D:\\Dmelanogaster\\Phosphoproteome\\new_data1_correct.txt"));
			
			String s1;
			Map<String,Double> hm2yw3 = new HashMap<String,Double>();
			while((s1 = br2.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String id = sp1[0];
				String pos = sp1[1];
				String site = id + "\t" + pos;
				double yw3 = Double.parseDouble(sp1[2]);
				if(yw3 != 0.0){
					hm2yw3.put(site, yw3);
				}
				
			}
			
			while((s1 = br1.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String id = sp1[0];
				String pos = sp1[1];
				String site = id + "\t" + pos;
				if(!hm2yw3.containsKey(site)){
					continue;
				}else{
					double yw3 = Double.parseDouble(sp1[3]);
					if(yw3 == 0.0){
						continue;
					}
					double foldchange = hm2yw3.get(site)/yw3;
					bw1.write(site);
					for(int i=2;i<sp1.length;i++){
						if(i == 3){
							continue;
						}
						double inten = Double.parseDouble(sp1[i]);
						double newinten = inten * foldchange;
						bw1.write("\t" + String.format("%.2f", newinten));
					}
					bw1.newLine();
				}
			}
			
			bw1.flush();
			bw1.close();
		}*/
	}
}
