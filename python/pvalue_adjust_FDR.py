import warnings
def pvalue_adjust(pvalue:[],method='FDR'):
    if method!='FDR':
        warnings.warn('Not support other Method, try to adjust p value by FDR!')
    leng=len(pvalue)
    pvalue_idx=[(i,pvalue[i]) for i in range(leng)]
    sortpvalue=sorted(pvalue_idx,key=lambda x:x[1])
    bh_fdr=[]
    if sortpvalue[-1][1]>1:
        bh_fdr.append((sortpvalue[-1][0],1))
    else:
        bh_fdr.append(sortpvalue[-1])
    for i in range(2,leng+1):
        rank=leng - i+1
        pval_idx= sortpvalue[leng - i]
        fdr=pval_idx[1]* leng /rank
        fdr_front=bh_fdr[-1][1]
        if fdr>1:
            bh_fdr.append((pval_idx[0],1))
        elif fdr>fdr_front:
            bh_fdr.append((pval_idx[0],fdr_front))
        else:
            bh_fdr.append((pval_idx[0],fdr))
    return [x[1] for x in sorted(bh_fdr,key=lambda x:x[0])]

if __name__=='__main__':
    a=[7.25E-09,1.89E-08,3.51E-08,9.09E-08,0.000000112,0.000000174,0.000000187,0.000000227,0.000000277,0.000000483,0.000000816]
    ar=pvalue_adjust(a)
    for i in range(len(a)):
        print('%s, %s'%(a[i],ar[i]) )

