import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class PhosNormByPro {
	public static void main(String[] args) throws IOException{
		File f11 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data1_correct.txt");
		File f12 = new File("D:\\Dmelanogaster\\Proteome\\new_data1_correct.txt");
		File f21 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data2.txt");
		File f22 = new File("D:\\Dmelanogaster\\Proteome\\new_data2.txt");
		File f31 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data3.txt");
		File f32 = new File("D:\\Dmelanogaster\\Proteome\\new_data3.txt");
		File f41 = new File("D:\\Dmelanogaster\\Phosphoproteome\\new_data4.txt");
		File f42 = new File("D:\\Dmelanogaster\\Proteome\\new_data4.txt");
		
		if(f11.exists() && f12.exists()){
			BufferedReader br11 = new BufferedReader(new FileReader(f11));
			BufferedReader br12 = new BufferedReader(new FileReader(f12));
			BufferedReader br21 = new BufferedReader(new FileReader(f21));
			BufferedReader br22 = new BufferedReader(new FileReader(f22));
			BufferedReader br31 = new BufferedReader(new FileReader(f31));
			BufferedReader br32 = new BufferedReader(new FileReader(f32));
			BufferedReader br41 = new BufferedReader(new FileReader(f41));
			BufferedReader br42 = new BufferedReader(new FileReader(f42));
			
			
			BufferedWriter bw1 = new BufferedWriter(new FileWriter("D:\\Dmelanogaster\\Normed\\KAA_wt_16.txt"));
			BufferedWriter bw2 = new BufferedWriter(new FileWriter("D:\\Dmelanogaster\\Normed\\KAA_pe_16.txt"));
			
			Map<String,List<Double>> data1 = tools.norm(br11,br12,6);
			Map<String,List<Double>> data2 = tools.norm(br21,br22,10);
			Map<String,List<Double>> data3 = tools.norm(br31,br32,8);
			Map<String,List<Double>> data4 = tools.norm(br41,br42,8);
			
			bw1.write("\tyw0\tyw3\tyw6\tyw9\tyw12\tyw15\tyw18\tyw21\tyw24\tyw27\tyw30\tyw33\tyw36\tyw39\tyw42\tyw45");
			bw1.newLine();
			bw2.write("\tpe0\tpe3\tpe6\tpe9\tpe12\tpe15\tpe18\tpe21\tpe24\tpe27\tpe30\tpe33\tpe36\tpe39\tpe42\tpe45");
			bw2.newLine();
			Iterator<String> it1 = data1.keySet().iterator();
			while(it1.hasNext()){
				String inf = it1.next();
				System.out.println(inf);
				if(data2.containsKey(inf) && data3.containsKey(inf)){
					List<Double> l1 = data1.get(inf);
					List<Double> l2 = data2.get(inf);
					List<Double> l3 = data3.get(inf);
					if(l1.get(0)*l2.get(0)*l1.get(1)*l2.get(1)*l1.get(2)*l2.get(2)*l1.get(3)*l2.get(3)* 
								l3.get(0)*l3.get(1)*l3.get(2)*l3.get(3)* l3.get(4)*l3.get(5)*l3.get(6)*l3.get(7) != 0.0){
						bw1.write(inf + "\t" + l1.get(0) + "\t" + l2.get(0) + "\t" + l1.get(1) + "\t" + l2.get(1) + "\t" + 
								l1.get(2) + "\t" + l2.get(2) + "\t" + l1.get(3) + "\t" + l2.get(3) + "\t" + 
								l3.get(0) + "\t" + l3.get(1) + "\t" + l3.get(2) + "\t" + l3.get(3) + "\t" + 
								l3.get(4) + "\t" + l3.get(5) + "\t" + l3.get(6) + "\t" + l3.get(7));
						bw1.newLine();
					}
					
				}
				if(data2.containsKey(inf) && data4.containsKey(inf)){
					List<Double> l1 = data1.get(inf);
					List<Double> l2 = data2.get(inf);
					List<Double> l4 = data4.get(inf);
					if(l2.get(4)*l2.get(5)*l1.get(4)*l2.get(6)*l2.get(7)*l2.get(8)*l1.get(5)*l2.get(9)* 
							l4.get(0)*l4.get(1)*l4.get(2)*l4.get(3)*l4.get(4)*l4.get(5)*l4.get(6)*l4.get(7) != 0.0){
						bw2.write(inf + "\t" + l2.get(4) + "\t" + l2.get(5) + "\t" + l1.get(4) + "\t" + l2.get(6) + "\t" + 
								l2.get(7) + "\t" + l2.get(8) + "\t" + l1.get(5) + "\t" + l2.get(9) + "\t" + 
								l4.get(0) + "\t" + l4.get(1) + "\t" + l4.get(2) + "\t" + l4.get(3) + "\t" + 
								l4.get(4) + "\t" + l4.get(5) + "\t" + l4.get(6) + "\t" + l4.get(7));
						bw2.newLine();
					}
					
				}
			}
			
			bw1.flush();
			bw1.close();
			bw2.flush();
			bw2.close();
		}
	}
}
