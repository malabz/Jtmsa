'''
MSA---star alignment
author: Juntao Chen
'''

from PSA_Kband import *
import numpy as np
from Extract_data import read_fasta

def insertGap(mark, seq):
    res = ""
    length = len(mark)
    for i in range(length):
        res += "-" * mark[i]
        if i < length - 1: 
            res += seq[i]
    return res

def MSA_star(S):
    # 1. to find the center sequence
    s_psa = [[-float('Inf')]*len(S) for _ in range(len(S))]
    
    print("-----------RUN-----------")
    print("Loading", end = "")
    
    for i in range(len(S)):
        if i % (len(S)//10) == 0.0:
            print(" => ", end="", flush=True)
        
        for j in range(len(S)):
            if j > i:
                tmp, _, _ = PSA_AGP_Kband(S[i], S[j])
                s_psa[i][j] = s_psa[j][i] = tmp
            elif j == i:
                s_psa[i][i] = 0

    C = np.argmax(np.sum(s_psa, axis=0))

    print("Loaded")
    print("center seq:", ','.join(S[C]))

    # 2. do pairwise alignments
    Strings = []
    for i in range(len(S)):
        if i != C:
            _, tmp1, tmp2 = PSA_AGP_Kband(S[C], S[i])
            Strings.append([tmp1, tmp2])

    # pprint(Strings)
    # 3. build the multiple alignment
    S_aligned = []
    markInsertion = [0]*(len(S[C]) + 1)
    for str2 in Strings:
        i = 0
        counter = 0
        for c in str2[0]:
            if c == '-': counter += 1
            else:
                markInsertion[i] = max(markInsertion[i], counter)
                counter = 0
                i += 1
            markInsertion[i] = max(markInsertion[i], counter)
    S_aligned = [""]*(len(S))
    S_aligned[C] = insertGap(markInsertion, S[C])
    idx = 0
    for str2 in Strings:
        mark = [0]*(len(str2[0])+1)
        total = 0
        pi = 0
        pj = 0
        for c in str2[0]:
            if c == '-': 
                total += 1
            else:
                mark[pi] = markInsertion[pj] - total
                pi += 1
                pj += 1
                while total != 0:
                    pi += 1
                    total -= 1
        mark[pi] = markInsertion[pj] - total
        if idx >= C:
            S_aligned[idx + 1] = insertGap(mark, str2[1])
        else:
            S_aligned[idx] = insertGap(mark, str2[1])
        idx += 1

    # 4. compute the SP value
    Value_SP = 0
    for i in range(len(S)):
        for j in range(len(S)):
            if j > i:
                Value_SP += Compute_two(S_aligned[i], S_aligned[j])

    print("SP : ", Value_SP)
    for str in S_aligned:
        print(' '.join(str))
    
    print("-----------END-----------")

    return Value_SP, S_aligned

if __name__ == "__main__":
    strs = read_fasta("dna500.fasta")
    MSA_star(strs)
