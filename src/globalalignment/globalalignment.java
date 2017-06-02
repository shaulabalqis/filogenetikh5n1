package globalalignment;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.Scanner;
import java.lang.Math.*;
//import UI.GlobalAlignmentUI;
//import UI.hasil;

/**
 *
 * @author Kiky
 */

public class globalalignment {
    public static class fungsi{
	
    //public static GlobalAlignmentUI prim=new UI.GlobalAlignmentUI();
    public static int[][] pam250;
    public static int[][] blosum62;
    public static int [][] dataTest;
    public static String dataQ; public static String[] DataQ;
    public static String dataD; public static String[] DataD;
    public static String[] sinkron={"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};
    public static int gapInisialisasi,gapIterasi,baris, kolom, pilihangap, pilihanmatriks,max;
    public static enum Status{kiri, atas, diagonal};
    public static Status[][] status;
    public static String hasilAlignQ, hasilAlignD;
    public static void inisialisasi(){
        try{
           pam250=new int [20][20];
           blosum62=new int [20][20];
           FileInputStream in=new FileInputStream("pam250.txt");
           Scanner baca = new Scanner(in);
           BufferedReader line=new BufferedReader(new InputStreamReader(in));
           for (int i=0; i<20; i++){
               for (int j=0; j<20; j++)
                   pam250[i][j]=baca.nextInt();
               line.readLine();
           }
           FileInputStream in2=new FileInputStream("blosum62.txt");
           Scanner baca2 = new Scanner(in2);
           BufferedReader line2=new BufferedReader(new InputStreamReader(in2));
           for (int i=0; i<20; i++){
               for (int j=0; j<20; j++)
                   blosum62[i][j]=baca2.nextInt();
               line2.readLine();
           }
        }
        catch(Exception e){
            System.out.println(e);
        }
    }
    public static void matrikstest(){
        char[] dataQ2=dataQ.toCharArray();
        char[] dataD2=dataD.toCharArray();
        baris=dataD2.length+1;
        kolom=dataQ2.length+1;
        dataTest=new int [baris][kolom];
        for (int i=0; i<kolom; i++)
            dataTest[0][i]=i*(-1)*gapInisialisasi;
        for (int i=0; i<baris; i++)
            dataTest[i][0]=i*(-1)*gapInisialisasi;
        DataQ=new String [kolom]; DataD=new String[baris];
        String temp="";
        DataQ[0]=null;DataD[0]=null;
        int i=1;
        //System.out.print("\t");
        for (int k=0; k<dataQ2.length;k++){
		temp += dataQ2[k];
                DataQ[i]=temp;
		//System.out.print (DataQ[i]+"\t");
		temp="";
                i++;
	}
	//System.out.println();
	i=1;
        for (int k=0; k<dataD2.length;k++){
                temp += dataD2[k];
                DataD[i]=temp;
		//System.out.print (DataD[i]+"\n");
		temp="";
                i++;
        }
	//System.out.println();
    }
    public static void algoritmaPAM250(){
        status=new Status[baris][kolom];
        status[0][0]=null;
        for (int i=1; i<baris; i++){
            status[i][0]=null;
            for (int j=1; j<kolom; j++){
                status[0][j]=null;
                int q=0,d=0,r;
                for (int k=0; k<20; k++){
                    if (DataQ[j].equals(sinkron[k]))
                        q=k;
                    if (DataD[i].equals(sinkron[k]))
                        d=k;
                }
                //System.out.print(d+","+q+"\t");
                r=pam250[d][q];
                //System.out.print(r+"\t");
                if (((dataTest[i-1][j-1]+r)>(dataTest[i][j-1]-gapIterasi))&&((dataTest[i-1][j-1]+r)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i-1][j-1]+r;
                    status[i][j]=Status.diagonal;
                }
                else if(((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j-1]+r))&&((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i][j-1]-gapIterasi;
                    status[i][j]=fungsi.Status.kiri;
                }
                else{
                    dataTest[i][j]=dataTest[i-1][j]-gapIterasi;
                    status[i][j]=fungsi.Status.atas;
                }
            }
            //System.out.println();
        }
        /* for (int i=0; i<baris; i++){
            for (int j=0; j<kolom; j++){
               if(status[i][j]==fungsi.Status.diagonal)
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t");
                else
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t\t");
            }
            System.out.println();
         }*/   
    }
    public static void algoritmaBlosum62(){
        status=new Status[baris][kolom];
        for (int i=1; i<baris; i++){
            for (int j=1; j<kolom; j++){
                int q=0,d=0,r;
                for (int k=0; k<20; k++){
                    if (DataQ[j].equals(sinkron[k]))
                        q=k;
                    if (DataD[i].equals(sinkron[k]))
                        d=k;
                }
                //System.out.print(d+","+q+"\t");
                r=blosum62[d][q];
                //System.out.print(r+"\t");
                if (((dataTest[i-1][j-1]+r)>(dataTest[i][j-1]-gapIterasi))&&((dataTest[i-1][j-1]+r)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i-1][j-1]+r;
                    status[i][j]=Status.diagonal;
                }
                else if(((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j-1]+r))&&((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i][j-1]-gapIterasi;
                    status[i][j]=Status.kiri;
                }
                else{
                    dataTest[i][j]=dataTest[i-1][j]-gapIterasi;
                    status[i][j]=Status.atas;
                }
            }
            //System.out.println();
        }
         /*for (int i=0; i<baris; i++){
            for (int j=0; j<kolom; j++){
                if(status[i][j]==Status.diagonal)
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t");
                else
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t\t");
            }
            System.out.println();
         }*/   
    }
    public static void alignment(){
        int d=baris-1; int q=kolom-1;max=dataTest[d][q];
       for (int j=1; j<kolom; j++){
            if (dataTest[baris-1][j]>max){
                max=dataTest[baris-1][j];
                d=baris-1;
                q=j;
            }
        }
        for (int j=1; j<baris; j++){
            if (dataTest[j][kolom-1]>max){
                max=dataTest[j][kolom-1];
                d=j;
                q=kolom-1;
            }
        }
        max=dataTest[d][q];
        //System.out.println(max+" "+status[d][q]+" "+d+","+q);
        boolean stop=false;
        hasilAlignD="";hasilAlignQ="";String gap="-";String space="\t";
        if (q!=kolom-1){
            for (int i=kolom-1; i>q; i--){
                hasilAlignQ=DataQ[i].toUpperCase()+space+hasilAlignQ;
                hasilAlignD=gap+space+hasilAlignD;
            }
        }
        else if (d!=baris-1){
            for (int i=baris-1; i>d; i--){
                hasilAlignD=DataD[i].toUpperCase()+space+hasilAlignD;
                hasilAlignQ=gap+space+hasilAlignQ;
            }
        }
        while (stop==false){
            if (status[d][q]==Status.diagonal){
                hasilAlignQ=DataQ[q].toUpperCase()+space+hasilAlignQ;
                hasilAlignD=DataD[d].toUpperCase()+space+hasilAlignD;
                d=d-1;q=q-1;
            }
            else if (status[d][q]==Status.atas){
                hasilAlignQ=gap+space+hasilAlignQ;
                hasilAlignD=DataD[d].toUpperCase()+space+hasilAlignD;
                d=d-1;
            }
            else{
                hasilAlignQ=DataQ[q].toUpperCase()+space+hasilAlignQ;
                hasilAlignD=gap+space+hasilAlignD;
                q=q-1;
            }
            if(DataQ[q]==null){
                for (int i=d; i>0; i--){
                    hasilAlignQ=gap+space+hasilAlignQ;
                    hasilAlignD=DataD[i].toUpperCase()+space+hasilAlignD;
                }
                stop=true;
            }
            else if (DataD[d]==null){
                for (int i=q; i>0; i--){
                    hasilAlignQ=DataQ[i].toUpperCase()+space+hasilAlignQ;
                    hasilAlignD=gap+space+hasilAlignD;
                }
                stop=true;
            }
        }
	
        //System.out.println("q :\t"+hasilAlignQ+"\nd :\t"+hasilAlignD);
    }
}
    
}
