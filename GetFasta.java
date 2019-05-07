import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class GetFasta {
	public static void main(String[] args) throws IOException{
		File f1 = new File("D:\\Dmelanogaster\\normed\\id_all.txt");
		File f2 = new File("D:\\Dmelanogaster\\ref_database\\fly_uni_ref.fasta");
		if(f1.exists()){
			BufferedReader br1 = new BufferedReader(new FileReader(f1));
			BufferedReader br2 = new BufferedReader(new FileReader(f2));
			BufferedWriter bw1 = new BufferedWriter(new FileWriter
					("D:\\Dmelanogaster\\normed\\phospho_seq_all.fa"));
			String s1;
			Set<String> hs1 = new HashSet<String>();
			while((s1 = br1.readLine()) != null){
				String[] sp1 = s1.split("_");
				String upid = sp1[0];
				hs1.add(upid);
			}
			
			String s2;
			boolean bl = false;
			Set<String> hs2 = new HashSet<String>();
			while((s2 = br2.readLine()) != null){
				if(s2.startsWith(">")){
					String[] sp2 = s2.split("\\|");
					String upid = sp2[1];
					if(hs1.contains(upid)){
						hs2.add(upid);
						bl = true;
						bw1.write(">" + upid);
						bw1.newLine();
					}else{
						bl = false;
					}
				}else{
					if(bl == true){
						bw1.write(s2);
						bw1.newLine();
					}else{
						continue;
					}
				}
			}
			bw1.flush();
			bw1.close();
		}
	}
}
