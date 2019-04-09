## convert IMPUT2 (.gen/.map) format to PLINK .ped (.map/.ped) format
import os,gzip,time
GEN_dir='/home/zm/main/schistosome/SNParray/LMX_impute2'
Imput_sample='/home/zm/main/schistosome/SNParray/LMX_impute2/chr1/chr1_new.sample'
PLINK_dir='/public1/tmp/xc/python'
plink_sep=' '

def log(content):
    print(time.strftime("[%Y-%m-%d %H:%M:%S] ", time.localtime(time.time())) + content)
def convert(Gen_files):
    ## write to .MAP of PLINK
    if not os.path.exists(PLINK_dir):
        os.makedirs(PLINK_dir)
    map_bw = open(os.path.join(PLINK_dir, 'convert.map'), 'w')
    ind_geno=[]
    al_idx=3;a2_idx=4 ## allele index in .GEN file
    log('start to read .gen of IMPUT2 ...')
    for gen in Gen_files:
        chrid=os.path.basename(gen).split('.')[0].replace('chr','')
        if not chrid=='22':
            continue
        log('chr%s start ...'%(chrid))
        br=gzip.open(gen,'r')
        # br.readline()## skip header
        limit_c=0
        for line in br:
            limit_c+=1
            # if limit_c>limit:
            #     break
            arr=line.decode().strip().split(' ')
            snp=arr[1].split(':')[0]
            if snp.startswith('rs'):
                snpid=snp
            else:
                snpid=chrid+':'+arr[2]
            map_bw.write(plink_sep.join([chrid,snpid,'0',arr[2]])+'\n')
            map_bw.flush()
            ind_size=int((len(arr)-5)/3)
            if len(ind_geno)==0:
                # ind_geno=[[] for x in range((len(arr)-5)/3)] ## init the ind_geno
                ind_geno=[[] for x in range(ind_size)] ## init the ind_geno

            for idx in range(ind_size):
                # geno_p=[float(arr[i]) for i in range(idx*3+5,(idx+1)*3+5)]
                # genotype=geno_p.index(max(geno_p))
                maxV=-1
                maxID=-1
                j=-1
                for i in range(idx*3+5,(idx+1)*3+5):
                    j+=1
                    V=float(arr[i])
                    if V>maxV:
                        maxV=V
                        maxID=j
                if maxID==0:
                    plink_gp=arr[al_idx]+plink_sep+arr[al_idx]
                elif maxID==1:
                    plink_gp = arr[al_idx] + plink_sep+arr[a2_idx]
                else:
                    plink_gp=arr[a2_idx]+plink_sep+arr[a2_idx]
                # ind_geno[idx].append(plink_gp)
                ind_geno[idx].append(plink_gp)
        br.close()
    map_bw.close()
    log('finish writing to .map of PLINK ...')
    ## write to PED
    log('start to write to .ped of PLINK ...')
    ped_bw = open(os.path.join(PLINK_dir, 'convert.ped'), 'w')
    ## read IMPUT sample file
    br=open(Imput_sample,'r')
    head_row=2
    for c in range(head_row):## skip header
        br.readline()
    idx=-1
    for line in br:
        idx+=1
        arr=line.strip().split(' ')
        ped_bw.write(plink_sep.join([arr[0].replace('MX2_F','')]+[arr[i] for i in [1,3,4,5,6]])+plink_sep+plink_sep.join(ind_geno[idx])+'\n')
        ped_bw.flush()
    br.close()
    ped_bw.close()
    log('all finished ...')


def run():
    Gen_files=[]
    for file in os.listdir(GEN_dir):
        if not (file.startswith('chr') and os.path.isdir(os.path.join(GEN_dir,file))):
            continue
        gen=os.path.join(GEN_dir, file, file + '.gen.gz')
        if not os.path.exists(gen):
            print('no such file: %s'%(gen))
            continue
        Gen_files.append(gen)
    Gen_files.sort(key=lambda x:int(os.path.basename(x).split('.')[0].replace('chr','')))
    convert(Gen_files)

if __name__=='__main__':
    run()
