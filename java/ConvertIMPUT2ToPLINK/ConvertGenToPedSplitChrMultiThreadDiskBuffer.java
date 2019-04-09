/*
 * Copyright © 2018 SNPLIFE.COM
 * All Rights Reserved 
 * 
 */
package com.snplife.demo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.Channels;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CodingErrorAction;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

/**
 * convert IMPUT2 (.gen/.sample) format to PLINK (.map/.ped) format
 * @author xuechao
 */
public class ConvertGenToPedSplitChrMultiThreadDiskBuffer {
    String GEN_dir="/home/zm/main/schistosome/SNParray/LMX_impute2";
    String Imput_sample="/home/zm/main/schistosome/SNParray/LMX_impute2/chr1/chr1_new.sample";
    String PLINK_dir="/public1/tmp/xc/java_multi";
    String plink_sep=" ";
    SimpleDateFormat df = new SimpleDateFormat("[yyyy-MM-dd HH:mm:ss] ");
    int threadNum=10;
    String TMPdir="/public1/tmp/xc/tmp";
//    int limitc=0;
//    int limit=10000;
    public void log(String content){
        System.out.println(this.df.format(new Date())+content);
    }
    public String join(Object[]obj,String sep){
        StringBuilder strb = new StringBuilder();
        if(obj.length==0){
            return strb.toString();
        }
        for(int i=0;i<obj.length-1;i++){
            strb.append((String)obj[i]+sep);
        }
        strb.append((String)obj[obj.length-1]);
        return strb.toString();
    }
    
    public void convert(List<String> genPaths) throws IOException, Exception{
        Date startT=new Date();
        // write to .MAP of PLINK
        if(!new File(PLINK_dir).exists()){
            new File(PLINK_dir).mkdirs();
        }        
        if(!new File(TMPdir).exists()){
            new File(TMPdir).mkdirs();
        }
        Object bo;
        int nt=threadNum;
        ExecutorService executorService = Executors.newFixedThreadPool(nt);
        Set<Future> fs = new HashSet<Future>();
        for(String genPath:genPaths){
            while (fs.size() >= nt) {
                for (Future fu : fs) {
                    bo = new Object();
                    //waiting time !!!!
                    try {
                        bo = fu.get(1000, TimeUnit.MILLISECONDS);
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
                BufferedWriter mapBw;
                BufferedWriter pedBw;
                BufferedReader br;
                String line;
                String []arr;
                List<BufferedWriter> indiGenoBw;
                BufferedReader brTmp;
                @Override
                public void run() {
                    String chr = new File(genPath).getName().split("\\.")[0].replace("chr", "");
                    log(String.format("chr%s: start to read .gen of IMPUT2 ...", chr));
                    try {
                        br = getBufferedReader(genPath);
                        mapBw = new BufferedWriter(new FileWriter(new File(PLINK_dir + java.io.File.separator + String.format("chr%s.covert.map", chr))));
                        indiGenoBw = new ArrayList<BufferedWriter>();
                        while ((line = br.readLine()) != null) {
                            arr = line.split(" ", -1);
                            String snpID;
                            if (arr[1].split(":")[0].startsWith("rs")) {
                                snpID = arr[1].split(":")[0];
                            } else {
                                snpID = chr + ":" + arr[2];
                            }
                            mapBw.append(join(new String[]{chr, snpID, "0", arr[2]}, plink_sep) + "\n");
                            mapBw.flush();
                            int indSize = (arr.length - 5) / 3;
                            if (indiGenoBw.isEmpty()) {
                                for (int i = 0; i < indSize; i++) {
                                    indiGenoBw.add(new BufferedWriter(new FileWriter(new File(TMPdir + java.io.File.separator + String.format("chr%s.indi%s.ped.tmp", chr,i)))));
                                }
                            }
                            for (int i = 0; i < indSize; i++) {
                                int j = -1;
                                int maxID = -1;
                                float maxV = -1;
                                float V;
                                for (int k = i * 3 + 5; k < (i + 1) * 3 + 5; k++) {
                                    j++;
                                    V = Float.valueOf(arr[k]);
                                    if (V > maxV) {
                                        maxID = j;
                                        maxV = V;
                                    }
                                }
                                String pedGeno = null;
                                switch (maxID) {
                                    case 0:
                                        pedGeno = arr[3] + plink_sep + arr[3];
                                        break;
                                    case 1:
                                        pedGeno = arr[3] + plink_sep + arr[4];
                                        break;
                                    case 2:
                                        pedGeno = arr[4] + plink_sep + arr[4];
                                        break;
                                }
                                indiGenoBw.get(i).append(plink_sep + pedGeno);
                            }
                        }
                        for (int i = 0; i < indiGenoBw.size(); i++) {
                            indiGenoBw.get(i).close();
                        }
                        indiGenoBw=null;
                        br.close();
                        mapBw.close();
                        log(String.format("chr%s: finish writing .map!", chr));
                        pedBw = new BufferedWriter(new FileWriter(new File(PLINK_dir + java.io.File.separator + String.format("chr%s.covert.ped", chr))));
                        br = getBufferedReader(Imput_sample);
                        for (int i = 0; i < 2; i++) {
                            br.readLine();
                        }
                        int idx = -1;
                        while ((line = br.readLine()) != null) {
                            idx++;
                            arr = line.split(" ", -1);
                            brTmp=getBufferedReader(TMPdir + java.io.File.separator + String.format("chr%s.indi%s.ped.tmp", chr,idx));
                            pedBw.append(join(new String[]{arr[0].replace("MX2_F", ""), arr[1], arr[3], arr[4], arr[5], arr[6]}, plink_sep) + brTmp.readLine() + "\n");
                            pedBw.flush();
                            brTmp.close();
                            new File(TMPdir + java.io.File.separator + String.format("chr%s.indi%s.ped.tmp", chr,idx)).delete();
                        }
                        br.close();
                        pedBw.close();
                    } catch (Exception ex) {
                        log("exception!!!"+ex.getLocalizedMessage());
                        Logger.getLogger(ConvertGenToPedSplitChrMultiThreadDiskBuffer.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    log(String.format("chr%s: finish writing .ped!", chr));
                }
            }));
        }
        executorService.shutdown();
        executorService.awaitTermination(24, TimeUnit.DAYS);
        try{
            new File(TMPdir).delete();
        }catch(Exception e){
            log(String.format("fail to delete TMPdir: %s, because: %s",TMPdir,e.getLocalizedMessage()));
        }
        Date endT=new Date();
        long time=endT.getTime()-startT.getTime();
        log(String.format("All finished! (time: %s min %s sec)", time/(1000*60),(time/1000)%60));
    }
    
    public void run() throws Exception{
        List<String> genFiles=new ArrayList<String>();
        for(File chr:new File(GEN_dir).listFiles()){
            if(!(chr.isDirectory()&&chr.getName().startsWith("chr"))){
                continue;
            }
            String file=chr.getAbsolutePath()+java.io.File.separator+chr.getName()+".gen.gz";
            if(!new File(file).exists()){
                log(String.format("no such file: %s", file));
                continue;
            }
            genFiles.add(file);
        }
        genFiles.sort(new Comparator<String>(){
            public int num(String s){
                return Integer.valueOf(new File(s).getName().split("\\.")[0].replace("chr", ""));
            }
            public int compare(String o1, String o2) {
                return Integer.compare(num(o1),num(o2));
            }
        });
        convert(genFiles);
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
}
