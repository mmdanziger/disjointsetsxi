from scipy.stats import binom
from sys import argv
import numpy as np
try:
    import tqdm
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
    p,Qn =np.loadtxt(argv[1]).T
    Qp=np.zeros_like(Qn)

    with Pool(processes=cpu_count()) as pool:
        if hasbar:
            unordered_res =  list(tqdm.tqdm(pool.imap_unordered(Qp_reduced, p), total=len(Qp)))
            for idx,val in unordered_res:
                Qp[idx]=val
        else:
            res = pool.map_async(Qp_reduced,p)
            Qp = np.array(res.get())
        
    ofname = argv[1].rstrip(".txt") + "_conv.txt"
    np.savetxt(ofname, np.array(Qp))