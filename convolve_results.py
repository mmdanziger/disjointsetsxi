from scipy.stats import binom
from sys import argv
import numpy as np
from os.path import basename
from collections import Counter
try:
    import tqdm
    hasbar=True
except ImportError:
    hasbar=False

from multiprocessing import Pool,cpu_count

zero_thresh=1e-13

def Qp_reduced(p):
    global Qn,hasbar
    N=len(Qn)
    n = int(round(p*N))
    all_terms=[]
    conv_term = lambda n: binom.pmf(n,N,p)
    keep_going=True
    Qp_out=0
    Qp_out += conv_term(n)*Qn[n]
    i=0
    while keep_going:
        keep_going=False
        if n+i<N:
            fwd_term = conv_term(n+i)
            if fwd_term>zero_thresh:
                keep_going=True
                Qp_out+=fwd_term*Qn[n+i]
        if n-i>0:
            bwd_term = conv_term(n-i)
            if bwd_term > zero_thresh:
                keep_going=True
                Qp_out += bwd_term*Qn[n-i]
        i+=1
    
    return (n,Qp_out) if hasbar else Qp_out

if __name__ == "__main__":
    Qn_list_0 = []
    for fname in argv[1:]:
        try:
            pQn =np.loadtxt(fname)
            Qn_list_0.append(pQn)
            p,Qn = pQn.T
        except:
            print("%s wrong format"%fname)
    lengthcount = sorted(Counter(map(len,Qn_list_0)).items(),key=lambda x: x[1],reverse=True)
    
    print("Found files (filelength,no_of_files_of_length)")
    print(lengthcount)
    
    if len(lengthcount)>1:
        Qn_list = list(filter(lambda x: len(x) == lengthcount[0][0],Qn_list_0))
        print("Tossing %i files for insufficient data count"%(len(Qn_list_0)-len(Qn_list)))
    else:
        Qn_list = Qn_list_0
    
    Qn_list = np.array(Qn_list)
    p=Qn_list[:,:,0].mean(axis=0)  
    Qn=Qn_list[:,:,1].mean(axis=0)  
    Qnstd = Qn_list[:,:,1].std(axis=0) 
    Qp=np.zeros_like(Qn)
    
    with Pool(processes=cpu_count()) as pool:
        if hasbar:
            unordered_res =  list(tqdm.tqdm(pool.imap_unordered(Qp_reduced, p), total=len(Qp)))
            for idx,val in unordered_res:
                Qp[idx]=val
        else:
            res = pool.map_async(Qp_reduced,p)
            Qp = np.array(res.get())
        
    ofname = "conv_" + basename(argv[1])
    oarray = np.array([p,Qn,Qnstd,Qp])
    np.savetxt(ofname,oarray)
   
