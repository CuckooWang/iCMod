import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class CalCricadianSiteKa {
	public static void main(String[] args) throws IOException{
		File f5 = new File("D:\\Dmelanogaster\\normed\\gps_out.txt");
		File f6 = new File("D:\\Dmelanogaster\\ref_database\\Drosophila_melanogaster.BDGP6.pep.all.fa");
		
		File f0 = new File("D:\\Dmelanogaster\\normed\\id.txt");
		
		if(/*f1.exists() && f2.exists() && f3.exists() && f4.exists()*/ f5.exists()){
			BufferedReader br5 = new BufferedReader(new FileReader(f5));
			BufferedReader br6 = new BufferedReader(new FileReader(f6));
			BufferedReader br0 = new BufferedReader(new FileReader(f0));
			
			BufferedWriter bw1 = new BufferedWriter(new FileWriter
					("D:\\Dmelanogaster\\normed\\KA.txt"));
			
			String s1;
			
			Map<String,String> GN = new HashMap<String,String>();
			while((s1 = br6.readLine()) != null){
				if(s1.startsWith(">")){
					String[] sp1 = s1.split("\\s+");
					String id = sp1[0].substring(1);
					String gn = null;
					for(int i=1;i<sp1.length;i++){
						if(sp1[i].startsWith("gene_symbol")){
							gn = sp1[i].substring(12);
							GN.put(id, gn);
						}
					}
				}else{
					continue;
				}
			}
			
			Set<String> site = new HashSet<String>();
			Map<String,Set<String>> kinase = new HashMap<String,Set<String>>();
			while((s1 = br0.readLine()) != null){
				String[] sp1 = s1.split("_");
				String id = sp1[0] + "\t" + sp1[1];
				site.add(id);
			}
			
			kinase = readinf(kinase,site,br5);
			
			Iterator<String> it1 = kinase.keySet().iterator();
			while(it1.hasNext()){
				String ka = it1.next();
				String[] sp1 = ka.split("\t");
				String kaname = sp1[0];
				int num = kinase.get(ka).size();
				bw1.write(kaname + "_" + GN.get(kaname) + "\t" + num);
				bw1.newLine();
			}
			
			bw1.flush();
			bw1.close();
		}
	}
	
	public static Map<String,Set<String>> readinf(Map<String,Set<String>> kinase, Set<String> site, BufferedReader br1) throws IOException{
		String s1;
		while((s1 = br1.readLine()) != null){
			String[] sp1 = s1.split("\t");
			String id = sp1[0] + "\t" + sp1[1];
			String ka = sp1[4] + "\t" + sp1[3];
			if(site.contains(id)){
				if(!kinase.containsKey(ka)){
					Set<String> hs1 = new HashSet<String>();
					hs1.add(id);
					kinase.put(ka, hs1);
				}else{
					Set<String> hs1 = kinase.get(ka);
					hs1.add(id);
					kinase.put(ka, hs1);
				}
			}else{
				continue;
			}
		}
		
		return kinase;
	}
}
