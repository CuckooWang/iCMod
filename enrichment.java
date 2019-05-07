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

public class enrichment {
	public static void main(String[] args) throws IOException{
		String sp = "wt";
		String SP = "WT";
		File f0 = new File("D:\\Dmelanogaster\\normed\\KAA_" + sp + "_16.txt");
		File f1 = new File("D:\\Dmelanogaster\\normed\\" + SP + "_out_16\\KA.txt");
		File f2 = new File("D:\\Dmelanogaster\\normed\\gps_out.txt");
		File f3 = new File("D:\\Dmelanogaster\\normed\\" + SP + "_out_16\\id.txt");
		
		File f6 = new File("D:\\Dmelanogaster\\ref_database\\Drosophila_melanogaster.BDGP6.pep.all.fa");
		
		if(f1.exists()){
			BufferedReader br0 = new BufferedReader(new FileReader(f0));
			BufferedReader br1 = new BufferedReader(new FileReader(f1));
			BufferedReader br2 = new BufferedReader(new FileReader(f2));
			BufferedReader br3 = new BufferedReader(new FileReader(f3));
			BufferedReader br6 = new BufferedReader(new FileReader(f6));
			BufferedWriter bw1 = new BufferedWriter(new FileWriter
					("D:\\Dmelanogaster\\normed\\" + SP + "_out_16\\" + sp + "_enrichout.txt"));
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
			
			double M = 0.0;
			Set<String> allsite = new HashSet<String>();
			while((s1 = br0.readLine()) != null){
				String[] sp1 = s1.split("\t");
				if(sp1[0].equals("")){
					continue;
				}
				String[] sp2 = sp1[0].split("_");
				String site = sp2[0] + "\t" + sp2[1];
				if(!allsite.contains(site)){
					M ++;
					allsite.add(site);
				}
				
			}
			
			Set<String> cirsites = new HashSet<String>();
			while((s1 = br3.readLine()) != null){
				cirsites.add(s1.trim());
			}
			double N = 1.0 * cirsites.size();
			
			Map<String,Double> ns = new HashMap<String,Double>();
			while((s1 = br1.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String[] sp2 = sp1[0].split("_");
				String id = sp2[0];  
				double n = Double.parseDouble(sp1[1]);
				ns.put(id, n);
			}
			
			
			
			Map<String,Set<String>> ms = new HashMap<String,Set<String>>();
			while((s1 = br2.readLine()) != null){
				String[] sp1 = s1.split("\t");
				String site = sp1[0] + "\t" + sp1[1];
				if(!allsite.contains(site)){
					continue;
				}
				String ka = sp1[4];
				if(ms.containsKey(ka)){
					Set<String> hs1 = ms.get(ka);
					hs1.add(site);
					ms.put(ka, hs1);
				}else{
					Set<String> hs1 = new HashSet<String>();
					hs1.add(site);
					ms.put(ka, hs1);
				}
			}
			
			Iterator<String> it1 = ns.keySet().iterator();
			while(it1.hasNext()){
				String ka = it1.next();
				System.out.println(ka);
				double n = ns.get(ka);
				double m = ms.get(ka).size() * 1.0;
				
				double enri = (n/N)/(m/M);
				double pVal = 0.0;
				if(enri >= 1.0){
					pVal = tools.pValuecal(M, m, N, n);
				}else{
					pVal = tools.pValuecal2(M, m, N, n);
				}
				bw1.write(ka + "_" + GN.get(ka) + "\t" + enri + "\t" + pVal);
				bw1.newLine();
			}
			
			bw1.flush();
			bw1.close();
		}
	} 
	
}
