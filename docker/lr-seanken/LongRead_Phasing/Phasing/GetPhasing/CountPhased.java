package GetPhasing;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class CountPhased
{
public static void main(String args[]) throws Exception
{
ArrayList<String[]> ret=new ArrayList<String[]>();
File bamFile = new File(args[0]); 
File vcfFile=new File(args[1]);
int ind=Integer.parseInt(args[2]);
String saveFile=args[3];

FileWriter savFil = new FileWriter(saveFile); 

print("Create VCF HashMap");

HashMap<String, Boolean> map= HashMapVCF(vcfFile,ind);
System.out.println(map.size());

print("Read in long read data!");
SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
SAMRecordIterator r = sr.iterator();
int iter=0;
int Good=0;
int Bad=0;


HashMap<String, String> UMIMap=new HashMap<String, String>(); 


while(r.hasNext()) {
iter=iter+1;

if(iter%10000==0)
{
System.out.println(Good);
System.out.println(Bad);
System.out.println((float)Good/(float)(Good+Bad));
System.out.println(iter);
}

SAMRecord read=r.next();

String seq=read.getReadString();
String cbc=read.getStringAttribute("CB");
String umi=read.getStringAttribute("UB");
String gene=read.getStringAttribute("gn");
String[] geneList=gene.split(",");

gene=geneList[0];



//getReferencePositionAtReadPosition
int start=read.getStart();
String chrom=read.getContig();
//chrom="chr"+chrom;
int readLen=read.getReadLength();
int All1=0;
int All2=0;


for(int i=0;i<readLen;i=i+1)
{

int pos=read.getReferencePositionAtReadPosition(i+1);
pos=pos;

if(pos<1){continue;}


String base=Character.toString(seq.charAt(i));


String key=chrom+"_"+Integer.toString(pos)+"_"+base;


if(map.containsKey(key))
{
Boolean val=map.get(key);
if(val)
{
All1=All1+1;
}
else
{
All2=All2+1;
}

}


}

int Tot=All1+All2;
if(Tot<1){continue;}
float rat=(float)Math.max(All1,All2)/(float)Tot;

boolean Ambig=false;
if(rat>.95)
{
Good=Good+1;
}
else
{
Bad=Bad+1;
Ambig=true;
}



String res=umi+" "+cbc+" "+gene;
String val="None";
if(All1>All2)
{
val="All1";
}
if(All2>All1)
{
val="All2";
}

if(Ambig)
{
val="Ambig";
}

if(UMIMap.containsKey(res))
{
String cur=UMIMap.get(res);
if(!cur.equals(val))
{
UMIMap.put(res,"Ambig");
}
}
else
{
UMIMap.put(res,val);
}




}

System.out.println("Get Allele Counts");

HashMap<String,Integer> counts=new HashMap<String,Integer>();
Iterator<String> it=UMIMap.keySet().iterator();


while (it.hasNext()) {
String key=it.next();
String val=UMIMap.get(key);
String[] split=key.split(" ");
String cbc=split[1];
String gene=split[2];
//String allele=split[3];
String res=cbc+" "+gene+" "+val;
//savFil.write(res+ System.lineSeparator());
if(!counts.containsKey(res))
{
counts.put(res,0);
}

counts.put(res,counts.get(res)+1);

}

System.out.println("Save!");


Iterator<String> it_cnt=counts.keySet().iterator();

while (it_cnt.hasNext()) {
String key=it_cnt.next();
Integer val=counts.get(key);

String res=key+" "+Integer.toString(val);
savFil.write(res+ System.lineSeparator());

}

r.close();
sr.close();
savFil.close();

}


public static HashMap<String,Boolean> HashMapVCF(File vcfFile,int ind) throws Exception
{

Scanner sc = new Scanner(vcfFile); 

HashMap<String,Boolean> vcfMap=new HashMap<String,Boolean>();

int iter=0;

while (sc.hasNextLine()) 
{
iter=iter+1;

if(iter %100000==0)
{
System.out.println(iter);
}
String line=sc.nextLine();
if('#'==line.charAt(0))
{
continue;
}

String[] split_line=line.split("\\s+");


String chrom=split_line[0];
String pos=split_line[1];
String ref=split_line[3];
String alt=split_line[4];
int genoPos=9+ind;
String geno=split_line[genoPos];

String key1=chrom+"_"+pos+"_"+ref;
String key2=chrom+"_"+pos+"_"+alt;

Boolean val1=true;
Boolean val2=false;

if(geno.equals("0|1"))
{
vcfMap.put(key1,val1);
vcfMap.put(key2,val2);
}


if(geno.equals("1|0"))
{
vcfMap.put(key1,val2);
vcfMap.put(key2,val1);
}


}

return(vcfMap);

}


public static void print(String toPrint)
{

System.out.println(toPrint);

}

}

