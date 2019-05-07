import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class tools {
	public static void putinMap(Map<String,Set<String>> hm1, String time,String gene){
		if(!hm1.containsKey(time)){
			Set<String> hs1 = new HashSet<String>();
			hs1.add(gene);
			hm1.put(time, hs1);
		}else{
			Set<String> hs1 = hm1.get(time);
			hs1.add(gene);
			hm1.put(time, hs1);
		}
	}
	
	public static Map<String, List<String>> getFastaSeq(File f1) throws IOException{
		BufferedReader br1 = new BufferedReader(new FileReader(f1));
		
		String s1;
		Map<String,List<String>> seqs = new HashMap<String,List<String>>();
		String upid = "";
		List<String> seq = new ArrayList<String>();
		while((s1 = br1.readLine()) != null){
			if(s1.startsWith(">")){
				if(upid != ""){
					seqs.put(upid, seq);
					seq = new ArrayList<String>();
				}
				String[] sp0 = s1.split("\\|");
				upid = sp0[1];
				seq.add(s1.trim());
			}else{
				seq.add(s1.trim());
			}
		}
		
		seqs.put(upid, seq);
		return seqs;
	}
	
	public static Map<String,List<Double>> norm(BufferedReader br1,BufferedReader br2,int num) throws IOException{
		String s1;
		Map<String,List<Double>> pro = new HashMap<String,List<Double>>();
		while((s1 = br2.readLine()) != null){
			String[] sp1 = s1.split("\t");
			String id = sp1[0];
			List<Double> l = new ArrayList<Double>();
			for(int i=1; i<num+1;i++){
				l.add(Double.parseDouble(sp1[i]));
			}
			pro.put(id, l);
		}
		
		Map<String,List<Double>> site = new HashMap<String,List<Double>>();
		while((s1 = br1.readLine()) != null){
			String[] sp1 = s1.split("\t");
			String id = sp1[0];
			String index = sp1[1];
			String inf = id + "_" + index;
			if(!pro.containsKey(id)){
				continue;
			}else{
				List<Double> pr = pro.get(id);
				List<Double> l = new ArrayList<Double>();
				for(int i=2; i<num+2;i++){
					double oldint = Double.parseDouble(sp1[i]);
					double newint = oldint/pr.get(i-2);
					l.add(newint);
				}
				site.put(inf, l);
			}
		}
		
		return site;
	}
	
	public static Map<String,String> readFasta(File f1) throws IOException{
		BufferedReader br1 = new BufferedReader(new FileReader(f1));
		
		String s1;
		Map<String,String> gn = new HashMap<String,String>();
		Set<String> gns = new HashSet<String>();
		while((s1 = br1.readLine()) != null){
			if(s1.startsWith(">")){
				String[] sp0 = s1.split("\\|");
				String upid = sp0[1];
				String[] sp1 = s1.split("\\s+");
				String genename = null;
				for(int i=0;i<sp1.length;i++){
					if(sp1[i].startsWith("GN=")){
						String[] sp2 = sp1[i].split("=");
						String[] sp3 = sp2[1].split("\\\\");
						for(int j=0;j<sp3.length;j++){
							genename = sp3[j];
							gns.add(genename);
							if(genename.endsWith("-RA") || genename.endsWith("-RB") || genename.endsWith("-RC") || genename.endsWith("-RD")
									|| genename.endsWith("-RE") || genename.endsWith("-RF")){
								String temgenename = genename.substring(0, genename.length()-3);
								if(gns.contains(temgenename)){
									continue;
								}else{
									genename = temgenename;
									gns.add(genename);
								}
								
							}
							gn.put(genename, upid);
						}

					}
				}
				if(genename == null){
					System.out.println(upid + "with no gene name!");
				}
			}
		}
		
		return gn;
	}
	
	public static double pValuecal(double N, double n, double M, double m){
		tools ke = new tools();
		double pv = 0.0;
		for(double i=m;i<=n&&i<M;i++){
			pv += Math.exp((ke.combnumber(M, i) + ke.combnumber(N-M, n-i)) - ke.combnumber(N, n));
		}
		return pv;
	}
	
	public static double pValuecal2(double N, double n, double M, double m){
		tools ke = new tools();
		double pv = 0.0;
		for(double i=0;i<=m;i++){
			pv += Math.exp((ke.combnumber(M, i) + ke.combnumber(N-M, n-i)) - ke.combnumber(N, n));
		}
		return pv;
	}
	
	public double combnumber(double N, double n){
		tools ke = new tools();
		double comnum = ke.factorial(N) - (ke.factorial(n) + ke.factorial(N-n));
		return comnum;
	}
	
	public double factorial(double i){
		double sum = Math.log(1.0);
		for(double j=1;j<=i;j++){
			sum += Math.log(j);
		}
		return sum;
	}
}
