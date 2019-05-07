import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FPKM1 {
	public static void main(String[] args) throws IOException{
		File f1 = new File("D:\\Dmelanogaster\\Transcriptome\\pe1_tot\\gene_exp.diff");
		File f2 = new File("D:\\Dmelanogaster\\Transcriptome\\pe2_tot\\gene_exp.diff");
		File f3 = new File("D:\\Dmelanogaster\\Transcriptome\\yw1_tot\\gene_exp.diff");
		File f4 = new File("D:\\Dmelanogaster\\Transcriptome\\yw2_tot\\gene_exp.diff");
		File f0 = new File("D:\\Dmelanogaster\\ref_database\\fly_uni_ref.fasta");
		
		BufferedReader br0 = new BufferedReader(new FileReader(f0));
		BufferedReader br1 = new BufferedReader(new FileReader(f1));
		BufferedReader br2 = new BufferedReader(new FileReader(f2));
		BufferedReader br3 = new BufferedReader(new FileReader(f3));
		BufferedReader br4 = new BufferedReader(new FileReader(f4));
		BufferedWriter bw1 = new BufferedWriter(new FileWriter
				("D:\\Dmelanogaster\\Transcriptome\\fpkm1_database\\fpkm1_4.fasta"));
		String s1;
		
		Map<String,Set<String>> hm1 = new HashMap<String,Set<String>>();
		String prefix = "pe1";
		while((s1 = br1.readLine()) != null){
			if(s1.startsWith("test_id")){
				continue;
			}
			String[] sp1 = s1.split("\t");
			String[] genes = sp1[2].split(",");
			String time1 = sp1[4];
			String time2 = sp1[5];
			double fpkm1 = Double.parseDouble(sp1[7]);
			double fpkm2 = Double.parseDouble(sp1[8]);
			for(int i=0;i<genes.length;i++){
				String gene = genes[i];
				if(fpkm1 >= 1.0){
					tools.putinMap(hm1,prefix+time1,gene);
				}
				if(fpkm2 >= 1.0){
					tools.putinMap(hm1,prefix+time2,gene);
				}
			}
		}
		
		prefix = "pe2";
		while((s1 = br2.readLine()) != null){
			if(s1.startsWith("test_id")){
				continue;
			}
			String[] sp1 = s1.split("\t");
			String[] genes = sp1[2].split(",");
			String time1 = sp1[4];
			String time2 = sp1[5];
			double fpkm1 = Double.parseDouble(sp1[7]);
			double fpkm2 = Double.parseDouble(sp1[8]);
			for(int i=0;i<genes.length;i++){
				String gene = genes[i];
				if(fpkm1 >= 1.0){
					tools.putinMap(hm1,prefix+time1,gene);
				}
				if(fpkm2 >= 1.0){
					tools.putinMap(hm1,prefix+time2,gene);
				}
			}
		}
		
		prefix = "yw1";
		while((s1 = br3.readLine()) != null){
			if(s1.startsWith("test_id")){
				continue;
			}
			String[] sp1 = s1.split("\t");
			String[] genes = sp1[2].split(",");
			String time1 = sp1[4];
			String time2 = sp1[5];
			double fpkm1 = Double.parseDouble(sp1[7]);
			double fpkm2 = Double.parseDouble(sp1[8]);
			for(int i=0;i<genes.length;i++){
				String gene = genes[i];
				if(fpkm1 >= 1.0){
					Z_tools.putinMap(hm1,prefix+time1,gene);
				}
				if(fpkm2 >= 1.0){
					Z_tools.putinMap(hm1,prefix+time2,gene);
				}
			}
		}
		
		prefix = "yw2";
		while((s1 = br4.readLine()) != null){
			if(s1.startsWith("test_id")){
				continue;
			}
			String[] sp1 = s1.split("\t");
			String[] genes = sp1[2].split(",");
			String time1 = sp1[4];
			String time2 = sp1[5];
			double fpkm1 = Double.parseDouble(sp1[7]);
			double fpkm2 = Double.parseDouble(sp1[8]);
			for(int i=0;i<genes.length;i++){
				String gene = genes[i];
				if(fpkm1 >= 1.0){
					Z_tools.putinMap(hm1,prefix+time1,gene);
				}
				if(fpkm2 >= 1.0){
					Z_tools.putinMap(hm1,prefix+time2,gene);
				}
			}
		}
		
		Set<String> gene1 = hm1.get("pe2h0");
		Set<String> gene2 = hm1.get("pe2h3");
		Set<String> gene3 = hm1.get("pe2h6");
		Set<String> gene4 = hm1.get("pe2h9");
		Set<String> gene5 = hm1.get("pe2h12");
		Set<String> gene6 = hm1.get("pe2h15");
		Set<String> gene7 = hm1.get("pe2h18");
		Set<String> gene8 = hm1.get("pe2h21");
		//Set<String> gene9 = hm1.get("pe1h15");
		//Set<String> gene10 = hm1.get("pe1h21");
		
		Set<String> allgene = new HashSet<String>();
		allgene.addAll(gene1);
		allgene.addAll(gene2);
		allgene.addAll(gene3);
		allgene.addAll(gene4);
		allgene.addAll(gene5);
		allgene.addAll(gene6);
		allgene.addAll(gene7);
		allgene.addAll(gene8);
		/*allgene.addAll(gene9);
		allgene.addAll(gene10);*/
		
		Map<String,List<String>> seqs = tools.getFastaSeq(f0); 
		Map<String,String> gn = tools.readFasta(f0);
		Set<String> same = new HashSet<String>();
		Iterator<String> it1 = allgene.iterator();
		while(it1.hasNext()){
			String genename = it1.next();
			if(!gn.containsKey(genename)){
				continue;
			}else{
				String upid = gn.get(genename);
				if(!same.contains(upid)){
					same.add(upid);
					List<String> seq = seqs.get(upid);
					for(int i=0;i<seq.size();i++){
						bw1.write(seq.get(i));
						bw1.newLine();
					}
				}
			}
			
		}
		bw1.flush();
		bw1.close();
	}
}
