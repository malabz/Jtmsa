'''
multiple sequence alignment
- center star method

author: Juntao Chen
date: 12.29.2021
'''

import numpy as np
import time

from PSA_Kband_memorysaving import PSA_AGP_Kband
from FASTA import readfasta
from score import spscore


def insertGap(mark, seq):
    res = ""
    length = len(mark)
    for i in range(length):
        res += "-" * mark[i]
        if i < length - 1: 
            res += seq[i]
    return res

def findCenterSeq(strs:list):
    """
    to find the center sequence

    Returns:
        index of center sequence
    """
    s_psa = [[-float('Inf')]*len(strs) for _ in range(len(strs))]
    
    for i in range(len(strs)):
        if i % (len(strs)//10) == 0.0:
            print(" => ", end="", flush=True)
        
        for j in range(len(strs)):
            if j > i:
                tmp, _, _ = PSA_AGP_Kband(strs[i], strs[j])
                s_psa[i][j] = s_psa[j][i] = tmp
            elif j == i:
                s_psa[i][i] = 0

    idxC = np.argmax(np.sum(s_psa, axis=0))

    return idxC

def psa(strs:list, idxC:int):
    """
    align center sequence with others
    """
    strsAligned = []
    for i in range(len(strs)):
        if i != idxC:
            _, tmp1, tmp2 = PSA_AGP_Kband(strs[idxC], strs[i])
            strsAligned.append([tmp1, tmp2])
    return strsAligned


def getGapsLoc(strsAligned:list, markInsertion:list, idxC:int):
    """
    compute the gaps location
    """
    for str in strsAligned:
        i = 0
        counter = 0
        for c in str[0]:
            if c == '-': 
                counter += 1
            else:
                markInsertion[i] = max(markInsertion[i], counter)
                counter = 0
                i += 1
            markInsertion[i] = max(markInsertion[i], counter)
    return markInsertion

def insertSeqsGap(strsAligned:list, markInsertion:list, strs:list, idxC):
    """
    insert gaps to all sequences
    """
    S_aligned = [""]*(len(strs))
    S_aligned[idxC] = insertGap(markInsertion, strs[idxC])
    idx = 0
    for str2 in strsAligned:
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
        if idx >= idxC:
            S_aligned[idx + 1] = insertGap(mark, str2[1])
        else:
            S_aligned[idx] = insertGap(mark, str2[1])
        idx += 1
    return S_aligned

def MSA_star(strs):
    sTime = time.time()
    # 1. to find the center sequence
    print("-----------RUN-----------")
    print("Loading", end = "")
    idxC = 0
    print("Loaded")
    print("center seq:", ''.join(strs[idxC]))

    # 2. do pairwise alignments
    strsAligned = psa(strs, idxC)

    # 3. build the multiple alignment
    markInsertion = [0]*(len(strs[idxC]) + 1)
    markInsertion = getGapsLoc(strsAligned, markInsertion, idxC)
    strsAligned = insertSeqsGap(strsAligned, markInsertion, strs, idxC)
    
    # 4. compute the SP value
    Value_SP = spscore(strsAligned)
    eTime = time.time()
    print("Run time : %.2f s"%(eTime - sTime))
    print("SP : ", Value_SP)
    for str in strsAligned:
        print(' '.join(str))
    print("-----------END-----------")


    return Value_SP, strsAligned

if __name__ == "__main__":
    labels, strs = readfasta("dna500.fasta")
    MSA_star(strs)


