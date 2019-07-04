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

'''
This script takes a sequence of files with microcanonical results, averages them to obtain Qn
and convolves them with the binomial distribution to obtain Qp.

Usage: 
    python3 convolve_results.py N data_files*
where N is the system size and data_files* is bash-expanded to files with raw measurements of Qn in two columns: p, Qn
It outputs four rows:
p - the percolation occupation probability (inherited from inputs)
Qn - the "microcanonical" average
Qnstd - the "microcanonical" standard deviation
Qp - the convolved "canonical" values (same sampling rate as Qn -- this script does not use the convolution to interpolate)

The p vector needs to be equally spaced, but does not need to sample every single value from 1 to N.
For discussion of the theory see Newman and Ziff, Phys. Rev. E 64, 016706  (2001).
If you use this implementation, please cite: Danziger, Gross and Buldyrev, arXiv:1902.03708 (2019)
'''



def Qp_reduced(pQn_idx):
    global p,Qn,hasbar,N,step_size
    n = int(round(p[pQn_idx]*N))
    
    all_terms=[]
    conv_term = lambda n: binom.pmf(n,N,p[pQn_idx])
    keep_going=True
    Qp_out=0
    Qp_out += conv_term(n)*Qn[pQn_idx]
    i=0
    while keep_going:
        keep_going=False
        if n+i*step_size<N:
            fwd_term = conv_term(n+i*step_size)
            if fwd_term>zero_thresh:
                keep_going=True
                Qp_out+=fwd_term*Qn[pQn_idx+i]
        if n-i*step_size>0:
            bwd_term = conv_term(n-i*step_size)
            if bwd_term > zero_thresh:
                keep_going=True
                Qp_out += bwd_term*Qn[pQn_idx-i]
        i+=1
    
    return (pQn_idx,Qp_out) if hasbar else Qp_out

if __name__ == "__main__":
    Qn_list_0 = []
    N = int(argv[1])
    for fname in argv[2:]:
        try:
            pQn =np.loadtxt(fname)
            p,Qn = pQn.T
            Qn_list_0.append(pQn)
            
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
    
    step_sizes = Counter(np.diff(np.round(p*N).astype(int)))
    if len(step_sizes) >1:
        raise ValueError("Uneven steps in data, unable to calculate convolution")
    step_size = list(step_sizes.keys())[0]
    
    with Pool(processes=cpu_count()) as pool:
        index_vector = np.array(list(range(len(p)))).astype(int)
        if hasbar:
            unordered_res =  list(tqdm.tqdm(pool.imap_unordered(Qp_reduced, index_vector), total=len(Qp)))
            for idx,val in unordered_res:
                Qp[idx]=val
        else:
            res = pool.map_async(Qp_reduced,index_vector)
            Qp = np.array(res.get())
        
    ofname = "conv_" + basename(argv[1])
    oarray = np.array([p,Qn,Qnstd,Qp])
    np.savetxt(ofname,oarray)
   
