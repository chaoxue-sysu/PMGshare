/*
 * Copyright © 2018 SNPLIFE.COM
 * All Rights Reserved 
 * 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.nio.channels.Channels;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CodingErrorAction;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

/**
 *
 * @author xuechao
 */
public class CoverageCalculationV3Multi {
    String bamDir="/public1/data/projects/prof_hc/HC_exome_2015/bam_hguo";
    String samtoolsPath="/sdb1/tools/wgss/samtools/bin/samtools";
    String depthResultPath = "/home/jh/projects/RUNER/Autism/ExonCoverage/ResultV3Multi.txt";
    String pedPath = "/home/jh/projects/RUNER/Autism/ExonCoverage/autism_lmx_moved.ped";
    String geneExonRange = "/home/jh/projects/RUNER/Autism/ExonCoverage/GeneExonRegion.txt";
    
    SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    
    int nt=2;
    int geneLimit = 100000;
    int sampleLimit=4;
    int MAXDepth=8;
    
    Set<String>chrSet;
    Set<String> sampleID;
    
    public static void main(String[] args) throws Exception {
        CoverageCalculationV3Multi cc=new CoverageCalculationV3Multi();
        cc.calculate();
    }

    public void calculate() throws Exception{
        System.out.println("[" + df.format(new Date()) + "] "  + " start to calculate ... ");
        BufferedReader br;
        String line;
        String arr[];
        br = getBufferedReader(geneExonRange);
        List<String> geneRange=new ArrayList<String>();
        Map<String,int[]> geneDepth=new HashMap<String,int[]>();
        int cou = 0;
        while ((line = br.readLine()) != null) {
            cou++;
            if (cou > geneLimit) {
                break;
            }
            arr=line.split(" ");
            geneDepth.put(arr[0], new int[4]);
            int exonLen=0;
            for(int i=1;i<arr.length;i++){
                geneRange.add(arr[0]+"%"+i+"%"+ arr[i]);
                String []pos=arr[i].split(":")[1].split("-");
                exonLen+=Integer.valueOf(pos[1])-Integer.valueOf(pos[0]);
            }   
            geneDepth.get(arr[0])[0]=exonLen;
            for (int i = 1; i < 4; i++) {
                geneDepth.get(arr[0])[i]=0;
            }
        }
        br.close();
        geneRange.sort(new Comparator<String>(){
            @Override
            public int compare(String o1, String o2) {
                String chr1=o1.split("%")[2].split(":")[0];
                String chr2=o2.split("%")[2].split(":")[0];
                int start1=0;
                try{
                start1=Integer.valueOf(o1.split("%")[2].split(":")[1].split("-")[0]);
                }catch(Exception e){System.out.println(o1);System.exit(0);}
                int start2=Integer.valueOf(o2.split("%")[2].split(":")[1].split("-")[0]);
                int end1=Integer.valueOf(o1.split("%")[2].split(":")[1].split("-")[1]);
                int end2=Integer.valueOf(o2.split("%")[2].split(":")[1].split("-")[1]);
                if (chr1.equals("X")) {
                    chr1 = "23";
                }
                if (chr1.equals("Y")) {
                    chr1 = "24";
                }
                if (chr1.equals("M")) {
                    chr1 = "25";
                }
                if (chr2.equals("X")) {
                    chr2 = "23";
                }
                if (chr2.equals("Y")) {
                    chr2 = "24";
                }
                if (chr2.equals("M")) {
                    chr2 = "25";
                }
                int ichr1 = Integer.valueOf(chr1);
                int ichr2 = Integer.valueOf(chr2);
                if (ichr1 > ichr2) {
                    return 1;
                } else if (ichr1 < ichr2) {
                    return -1;
                } else {
                    if (start1 > start2) {
                        return 1;
                    } else if (start1 < start2) {
                        return -1;
                    } else {
                        if (end1 > end2) {
                            return 1;
                        } else if (end1 < end2) {
                            return -1;
                        } else {
                            return 0;
                        }
                    }
                }
            }
        });

        System.out.println("[" + df.format(new Date()) + "] " +"finish sorting...");
        System.out.println("[" + df.format(new Date()) + "] " +"start to make chr index...");
        List<CHRIndex> soredChrIndex=new ArrayList<CHRIndex>();
        // init the list
        soredChrIndex.add(new CHRIndex(geneRange.get(0).split("%")[2],new HashSet<String>()));
        soredChrIndex.get(0).pos.add(geneRange.get(0));
        for (int i=1;i<geneRange.size();i++){
            CHRIndex former=soredChrIndex.get(soredChrIndex.size()-1);
            String formerChrP=former.posIndex;
            String formerChr =formerChrP.split(":")[0];
            String[] formerPos =formerChrP.split(":")[1].split("-");
            int formerStart=Integer.valueOf(formerPos[0]);
            int formerEnd=Integer.valueOf(formerPos[1]);
            
            String presentPosS=geneRange.get(i);
            String presentChrP=presentPosS.split("%")[2];
            String presentChr=presentChrP.split(":")[0];
            if(presentChr.equals(formerChr)){
                String[] presentPos=presentChrP.split(":")[1].split("-");
                int presentStart=Integer.valueOf(presentPos[0]);
                int presentEnd=Integer.valueOf(presentPos[1]);
                if(formerStart<=presentStart && formerEnd>=presentStart){
                    former.pos.add(presentPosS);
                    if(formerEnd<presentEnd){
                        former.posIndex=formerChr+":"+formerStart+"-"+presentEnd;
                    }
                    continue;
                }
            }
            soredChrIndex.add(new CHRIndex(geneRange.get(i).split("%")[2], new HashSet<String>()));
            soredChrIndex.get(soredChrIndex.size()-1).pos.add(geneRange.get(i));
        }
        // get sample id
        sampleID=new HashSet<String>();
        br= getBufferedReader(pedPath);
        br.readLine();
        int scou=0;
        while((line=br.readLine())!=null){
            String sampleid=line.split("\t",-1)[1];
            if (!new File(new File(bamDir).getAbsolutePath()+java.io.File.separator+sampleid+"_final.bam").exists()){
                continue;
            }
            scou++;
            if(scou>sampleLimit){
                break;
            }
            sampleID.add(sampleid); 
        }
        br.close();
        System.out.println("sample size: "+sampleID.size());
        Object bo;
        ExecutorService executorService = Executors.newFixedThreadPool(nt);
        Set<Future> fs = new HashSet<Future>();
        for(String sample:sampleID){
            // multi thread
             while (fs.size() >= nt) {
                for (Future fu : fs) {
                    bo = new Object();
                    //waiting time !!!!
                    try {
                        bo = fu.get(10, TimeUnit.MILLISECONDS);
                    } catch (Exception e) {
                        continue;
                    }
                    if (bo == null) {
                        fs.remove(fu);
                        break;
                    }
                }
            }
            fs.add(executorService.submit(new Runnable() {
                @Override
                public void run() {
                    while(true){
                        System.out.println("[" + df.format(new Date()) + "] " + sample + " start... ");
                        // get stream from samtools output
                        BufferedReader coverbr=null;
                        String line;
                        String sampleBAM = new File(bamDir).getAbsolutePath() + java.io.File.separator + sample + "_final.bam";
                        String calculateDepth = samtoolsPath + " depth " + sampleBAM;
                        Process process = null;
                        try {
                            process = Runtime.getRuntime().exec(calculateDepth);
                            coverbr = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        } catch (Exception ex) {
                            try {
                                coverbr.close();
                            } catch (IOException ex1) {
                                Logger.getLogger(CoverageCalculationMulti.class.getName()).log(Level.SEVERE, null, ex1);
                            }
                            Logger.getLogger(CoverageCalculationMulti.class.getName()).log(Level.SEVERE, null, ex);
                            break;
                        }

                        int geneIdx = 0;
                        try {
                            while ((line = coverbr.readLine()) != null) {
                                if (geneIdx > soredChrIndex.size() - 1) {
                                    break;
                                }
                                String[] arr = line.split("\t");
                                String[] posIndex = soredChrIndex.get(geneIdx).posIndex.split(":");
                                String chr = posIndex[0];
                                int posStart = Integer.valueOf(posIndex[1].split("-")[0]);
                                int posEnd = Integer.valueOf(posIndex[1].split("-")[1]);
                                // synonym of mitochondrial chromosome, M or MT
                                if (arr[0].equals("MT")){
                                    arr[0]="M";
                                }
                                if(chrCompare(arr[0],chr)>0){
                                    geneIdx += 1;
                                    continue;
                                }
                                if (chrCompare(arr[0],chr)<0) {
                                    continue;
                                }
                                int Spos = Integer.valueOf(arr[1]);
                                if (Spos >= posStart) {
                                    if (Spos < posEnd) {
                                        for (String exonPos : soredChrIndex.get(geneIdx).pos) {
                                            String geneN = exonPos.split("%")[0];
                                            String[] exonPosArr = exonPos.split("%")[2].split(":");
                                            int exonPosStart = Integer.valueOf(exonPosArr[1].split("-")[0]);
                                            int exonPosEnd = Integer.valueOf(exonPosArr[1].split("-")[1]);
                                            if (Spos >= exonPosStart && Spos < exonPosEnd) {
                                                geneDepth.get(geneN)[1] += 1;
                                                geneDepth.get(geneN)[3] += Integer.valueOf(arr[2]);
                                                if (Integer.valueOf(arr[2]) >= MAXDepth) {
                                                    geneDepth.get(geneN)[2] += 1;
                                                }
                                            }
                                        }
                                    } else {
                                        geneIdx += 1;
                                    }
                                }
                            }
                        } catch (IOException ex) {
                            Logger.getLogger(CoverageCalculationV3Multi.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        try {
                            process.waitFor();
                        } catch (InterruptedException ex) {
                            Logger.getLogger(CoverageCalculationMulti.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        process.destroy();
                        try {
                            coverbr.close();
                        } catch (IOException ex) {
                            Logger.getLogger(CoverageCalculationV3Multi.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        System.out.println("[" + df.format(new Date()) + "] " + sample + " finish... ");
                        break;
                    }
                }
            }));
        }
        executorService.shutdown();
        executorService.awaitTermination(24, TimeUnit.DAYS);
        System.out.println("[" + df.format(new Date()) + "] " + " start write to file... ");
        BufferedWriter bw=getBufferedWriter(depthResultPath,false);
        bw.write("GeneName\tExonLength\tExonCoverageMore_1\tExonCoverageMore_"+MAXDepth+"\tExonAvgDepth\n");
        for(String gene:geneDepth.keySet()){
            int count[]=geneDepth.get(gene);
            bw.write(gene+"\t"+count[0]+"\t"+Float.valueOf(count[1])/count[0]/sampleID.size()+"\t"+Float.valueOf(count[2])/count[0]/sampleID.size()+"\t"+Float.valueOf(count[3])/count[0]/sampleID.size()+"\n");
        }
        bw.close();
        System.out.println("[" + df.format(new Date()) + "] " + " all finish... ");
    }
    static private BufferedReader getBufferedReader(String filePath) throws Exception {
        BufferedReader br = null;
        File dataFile = new File(filePath);
        if (dataFile.exists()) {
            CharsetDecoder decoder = Charset.forName("UTF-8").newDecoder();
            decoder.onMalformedInput(CodingErrorAction.IGNORE);
            try {
                br = new BufferedReader(Channels.newReader(Channels.newChannel(new GZIPInputStream(Channels.newInputStream(new RandomAccessFile(dataFile, "r").getChannel()))), decoder, 1024 * 1024));
            } catch (ZipException ex) {
                br = new BufferedReader(Channels.newReader(Channels.newChannel(Channels.newInputStream(new RandomAccessFile(dataFile, "r").getChannel())), decoder, 1024 * 1024));
            }
        } else {
            throw new Exception("No input file: " + dataFile.getCanonicalPath());
        }
        return br;
    }
    static public BufferedWriter getBufferedWriter(String filePath, boolean isGzip) throws Exception {
        BufferedWriter bw = null;
        File dataFile = new File(filePath);
        if (!dataFile.getParentFile().exists()) {
            dataFile.getParentFile().mkdirs();
        }
        CharsetEncoder decoder = Charset.forName("UTF-8").newEncoder();
        decoder.onMalformedInput(CodingErrorAction.IGNORE);
        if (isGzip) {
            bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dataFile))));
        } else {
            bw = new BufferedWriter(new FileWriter(dataFile));
        }
        return bw;
    }

    private int chrCompare(String chr1, String chr2) {
        if (chr1.equals("X")) {
            chr1 = "23";
        }
        if (chr1.equals("Y")) {
            chr1 = "24";
        }
        if (chr1.equals("M")) {
            chr1 = "25";
        }
        if (chr2.equals("X")) {
            chr2 = "23";
        }
        if (chr2.equals("Y")) {
            chr2 = "24";
        }
        if (chr2.equals("M")) {
            chr2 = "25";
        }
        if (Integer.valueOf(chr1) > Integer.valueOf(chr2)) {
            return 1;
        } else if (Integer.valueOf(chr1) < Integer.valueOf(chr2)) {
            return -1;
        } else {
            return 0;
        }
    }
    class CHRIndex{
        String posIndex;
        Set<String>pos;
        public CHRIndex(String posIndex,Set<String>pos){
            this.pos=pos;
            this.posIndex=posIndex;
        }
    }
    
}
