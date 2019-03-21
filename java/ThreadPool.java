/*
 * Copyright © 2018 SNPLIFE.COM
 * All Rights Reserved 
 */
package com.snplife.test;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *  Thread pool demo code
 *
 * @author xuechao
 */
public class ThreadPool {
    SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    int threadNum=3;
    public static void main(String[] args) throws Exception {
        new ThreadPool().multiTask();
    }

    public void multiTask() throws Exception {
        Object bo;
        int nt=threadNum;
        ExecutorService executorService = Executors.newFixedThreadPool(nt);
        Set<Future> fs = new HashSet<Future>();
        for (int i=0;i<10;i++) {
            while(fs.size()>=nt){
                for(Future fu:fs){
                    bo = new Object();
                    //waiting time !!!!
                    try {bo = fu.get(1000, TimeUnit.MILLISECONDS);} catch(Exception e){continue;}
                    if (bo == null) {
                        fs.remove(fu);
                        break;
                    }
                }
            }
            fs.add(executorService.submit(new SingleTask(i)));
            System.out.println("["+df.format(new Date())+"] "+ i + " start... ");
        }
        executorService.shutdown();
        executorService.awaitTermination(24, TimeUnit.DAYS);
        System.out.println("all finish");
    }
    /**
    *  single task class
    *
    * @author xuechao
    */
    class SingleTask implements Runnable {
        int idx;
        public SingleTask(int i){
            this.idx=i;
        }

        @Override
        public void run() {
            String chrNum = "chr"+idx;
            Date startT = new Date();
            //System.out.println("["+df.format(startT)+"] "+ chrNum + " start... ");
            try {
                Thread.sleep((idx*2+1)*1000);
            } catch (InterruptedException ex) {
                Logger.getLogger(PopulationBasedVCFMulti.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            Date endT = new Date();
            System.out.println("["+df.format(endT)+"] "+ chrNum + " finish! (time: "+(endT.getTime()-startT.getTime())/(1000*60)+" min "+((endT.getTime()-startT.getTime())/1000)%60+" sec)");
        }
    }
}
