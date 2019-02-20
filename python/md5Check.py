#!/bin/python
import hashlib
import os, sys
def GetStrMd5(src):
    m0 = hashlib.md5()
    m0.update(src)
    print
    m0.hexdigest()
    pass
def GetBigFileMd5(filename):
    if not os.path.isfile(filename):
        return
    myhash = hashlib.md5()
    f = open(filename, 'rb')
    while True:
        b = f.read(8096)
        if not b:
            break
        myhash.update(b)
    f.close()
    return myhash.hexdigest()
def CalcSha1(filepath):
    with open(filepath, 'rb') as f:
        sha1obj = hashlib.sha1()
        sha1obj.update(f.read())
        hash = sha1obj.hexdigest()
        print(hash)
        return hash
def CalcMD5(filepath):
    with open(filepath, 'rb') as f:
        md5obj = hashlib.md5()
        md5obj.update(f.read())
        hash = md5obj.hexdigest()
        print(hash)
        return hash
def checkMD5ByPath(hashfile,md5):
    re=[]
    if not os.path.exists(hashfile):
        hashfile = os.path.join(os.path.dirname(__file__), hashfile)
        if not os.path.exists(hashfile):
            re=[md5,"404"]
            return re
        else:
            checked_md5=GetBigFileMd5(hashfile)
    else:
        checked_md5=GetBigFileMd5(hashfile)
    if checked_md5!=md5:
        re=[md5,checked_md5]
    else:
        re=[md5,"200"]
    return re
def getPathList(dirPath, md5Path):
    re=[]
    file=open(md5Path, "r")
    file.readline()
    for line in file:
        cells=line.strip().split("\t")
        re.append(([os.path.join(dirPath,cells[0],cells[1]),cells[2]],None))
    file.close()
    return re

def call_checkMD5(request,result):
    logDirPath="F:\Projects\TCGA\md5_checked"
    log=open(os.path.join(logDirPath,"all.log"),"a")
    error=open(os.path.join(logDirPath,"error.log"),"a")
    print(result)
    print(request.args[0])
    log.write("\t".join(result)+"\n")
    if not result[1] == "200":
        error.write("\t".join(result)+"\n")
    log.close()
    error.close()
def calcMD5MultiThread(dirPath,md5Path,logDirPath,threadNum,multiThread=True):
    if multiThread:
        poolsize = threadNum
        import threadpool as tp
        pool = tp.ThreadPool(poolsize)
        para=getPathList(dirPath,md5Path)
        requests = tp.makeRequests(checkMD5ByPath, para, call_checkMD5)
        [pool.putRequest(req) for req in requests]
        pool.wait()



if __name__ == "__main__":
    # if len(sys.argv) == 2:
    #     hashfile = sys.argv[1]
    # else:
    #     print("no filename")
    calcMD5MultiThread("F:\Projects\TCGA\WXS\LUAD",
                       "F:\Projects\TCGA\md5\LUAD_WXS_bam_gdc_manifest.2018-03-29.txt",
                       "F:\Projects\TCGA\md5_checked",40)
    calcMD5MultiThread("F:\Projects\TCGA\WXS\LUAD",
                       "F:\Projects\TCGA\md5\LUAD_WXS_bam_gdc_manifest.2018-03-29.txt",
                       "F:\Projects\TCGA\md5_checked",40)
