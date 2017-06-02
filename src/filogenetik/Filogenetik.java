package filogenetik;
/**
 *
 * @author Kiky
 */
import globalalignment.*;
import java.io.*;
import java.util.*;
//import javax.swing.JTable;
//import UI.FilogenetikUI;
public class Filogenetik {
    public static UI.FilogenetikUI prim =new UI.FilogenetikUI();
    public static int pilihanKlustering, jumlah;
    public static double[][] skor, skorBaru;//skorBar untuk setelah klustering
    public static String temp;
    public static String[] seq;
    public static double[] l ;//height
    public static int[] N, panjangSeq;
    public static String[] namaFile0={"ProteinFastaResultsKorea.fasta","ProteinFastaResultsThailand.fasta","ProteinFastaResultsGuangdong.fasta", "ProteinFastaResultsFrance.fasta","ProteinFastaResultsGermany.fasta"};
    public static String[] idFile0={"Korea", "Thai", "China", "France", "Germany"}, kluster0={"K", "T", "C", "F", "G"};
    public static String [] namaFile,idFile,kluster;
    public static Scanner in;
    public static BufferedReader baca;
    public static FileInputStream read;
    public static void init(){
	skor=new double[jumlah][jumlah]; 
	skorBaru=new double[jumlah][jumlah];
	seq=new String [jumlah];
	l=new double [jumlah] ;
	N=new int [jumlah];
	for (int i=0; i<jumlah; i++)
		N[i]=1;
	panjangSeq=new int[jumlah];
	namaFile=new String[jumlah];
	idFile=new String[jumlah];
	kluster=new String[jumlah];
    }	
    public static void bacaSequence(){
	try{
	    String tambah,file;
	    for (int i=0; i<jumlah; i++){
		    file=namaFile[i];
		temp="";
		read=new FileInputStream (file);
		baca=new BufferedReader(new InputStreamReader(read));
	        baca.readLine();
		tambah=baca.readLine();
	        while (tambah.length()>0){
			temp=temp+tambah;
			tambah=baca.readLine();
		}
		seq[i]=temp;
		panjangSeq[i]=temp.toCharArray().length;
	    }
	}
	catch(Exception e){
	    System.out.println(e);
	}
    }
    public static void alignment(){
	int k;
	for (int i=0; i<jumlah; i++){
	    for (int j=i+1; j<jumlah; j++){
		globalalignment.fungsi.dataQ=seq[i];globalalignment.fungsi.dataD=seq[j];
		globalalignment.fungsi.inisialisasi();
		globalalignment.fungsi.matrikstest();
                if (globalalignment.fungsi.pilihanmatriks==1){
                    globalalignment.fungsi.algoritmaPAM250();
                    globalalignment.fungsi.alignment();
                }
                else if (globalalignment.fungsi.pilihanmatriks==2){
                    globalalignment.fungsi.algoritmaBlosum62();
                    globalalignment.fungsi.alignment();
                }
		System.out.println("q("+idFile[i]+"):\t"+globalalignment.fungsi.hasilAlignQ+"\nd("+idFile[j]+"):\t"+globalalignment.fungsi.hasilAlignD+"\nskor: "+globalalignment.fungsi.max+"\n");
		skor[i][j]=globalalignment.fungsi.max;
		//skorBaru[i][j]=skor[i][j];
	    }
	}
        double temp=0;
        for (int i=0; i<jumlah; i++){
            for (int j=i+1; j<jumlah; j++){
                temp=temp+skor[i][j];
            }
        }
        for (int i=0; i<jumlah; i++){
            for (int j=i+1; j<jumlah; j++){
                skorBaru[i][j]=skor[i][j]/temp;
                skor[i][j]=skor[i][j]/temp;
            }
        }
    }
    public static void klusteringUPGMA(){
	double min=10000;int k=-1,m=-1;
	for (int i=0; i<jumlah; i++){
	    for (int j=i+1; j<jumlah; j++){
		if (skorBaru[i][j]<min){
		    min=skorBaru[i][j];
		    k=i;m=j;
		}		    
	    }
	}
	if (kluster[k]!=null&&kluster[m]!=null)
	    kluster[k]=kluster[k]+"-"+kluster[m];
	for (int i=0; i<jumlah; i++){
	    if (i!=k&&i!=m){
		if (k<i){
		    if (m<i)
			skorBaru[k][i]=(N[k]*skorBaru[k][i]+N[m]*skorBaru[m][i])/(N[k]+N[m]);
		    else if (i<m)
			skorBaru[k][i]=(N[k]*skorBaru[k][i]+N[m]*skorBaru[i][m])/(N[k]+N[m]);
		}
		else if (i<k){
		    if (m<i)
			skorBaru[k][i]=(N[k]*skorBaru[i][k]+N[m]*skorBaru[m][i])/(N[k]+N[m]);
		    else if (i<m)
			skorBaru[k][i]=(N[k]*skorBaru[i][k]+N[m]*skorBaru[i][m])/(N[k]+N[m]);
		}
	    }
	}
	N[k]=N[k]+N[m];
	kluster[m]=null;
	N[m]=0;
	l[k]=min/2;
	for (int j=0; j<jumlah; j++){
	    skorBaru[j][m]=Double.NaN;
	    skorBaru[m][j]=Double.NaN;
	}
	for (int i=0; i<jumlah; i++){
	    if (kluster[i]!=null)
		System.out.print("\t"+kluster[i]);
	    if (N[i]>1)
		System.out.print("("+l[i]+")");
	}
	System.out.print("\n"+kluster[0]);
	for (int j=1; j<jumlah; j++){
	    if (kluster[j]!=null)
	    System.out.print("\n"+kluster[j]);
	     for (int i=0; i<jumlah; i++){
		    if (i<j&&Double.isNaN(skorBaru[i][j])==false)
		    System.out.print("\t"+skorBaru[i][j]);
	     }
	}
	System.out.println("\n");    
    } 
    
    public static void cetakSistem(){
	System.out.print("DISTANCE :\n\t");
	for (int i=0; i<jumlah; i++){
	    System.out.print(idFile[i]+"\t");
	}
	System.out.print("\n"+idFile[0]);
	for (int j=1; j<jumlah; j++){
	    System.out.print("\n"+idFile[j]);
	     for (int i=0; i<jumlah; i++){
		    if (i<j)
		    System.out.print("\t"+skor[i][j]);
	     }
	}
	System.out.println("\n\nHASIL KLUSTER :\nNama Kluster(height)\n");
    
    }
    public static void main(String[] args) {
	prim.setVisible(true);
    }
}
