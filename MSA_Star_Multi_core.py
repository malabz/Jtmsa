'''
star alignment with multiprocessing
author: Juntao Chen
'''

import time
import numpy as np
import multiprocessing as mp

from PSA_Kband import PSA_AGP_Kband
from score import spscore
from FASTA import readfasta


def find_censeq_sub(name:str, param:list):
    """
    find the center sequence
    """
    strs = param[0]
    k = param[1]
    nums = param[2]
    length = len(strs)

    size = (length * (length - 1) // 2) // nums + 1
    res = [0] * size

    for i in range(length):
        for j in range(i+1, length):
            idx = getIdx(length, i, j)
            if (idx) % nums == k:
                score, _ = PSA_AGP_Kband(A = strs[i], B = strs[j], get_score = 1)
                res[idx//nums] = score

    return res

def getIdx(length:int, i:int, j:int):
    """
    Get the index of the location
    """
    return ((2*length - i - 1) * i // 2) + j - i


def insertGap(mark:list, seq:str):
    """
    insert gaps to sequences

    Args:
        mark: gaps loc list
        seq: sequence

    Returns:
        sequence
    """
    res = ""
    length = len(mark)
    for i in range(length):
        res += "-" * mark[i]
        if i < length - 1:
            res += seq[i]
    return res

def findCenterSeq(strs:list, num_cores:int):
    """
    find center sequence with multi cores
    """
    length = len(strs)
    s_psa = [[0]*length for _ in range(length)]
    param_dict = {}
    for i in range(num_cores):
        param_dict["task"+str(i)] = (strs, i, num_cores)
    pool = mp.Pool(num_cores)
    results = [pool.apply_async(find_censeq_sub, args=(name, param)) for name, param in param_dict.items()]
    results = [p.get() for p in results]
    print("=>>>>>", end='')
    for i in range(length):
        for j in range(i + 1, length):
            idx = getIdx(length, i, j)
            s_psa[i][j] = s_psa[j][i] = results[idx%num_cores][idx//num_cores]
    return np.argmax(np.sum(s_psa, axis=0))

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


def MSA_star_Multicores(strs, num_cores = -1):

    sTime = time.time()
    if num_cores == -1: num_cores = mp.cpu_count()
    if num_cores == None: raise ValueError("Can not get the number of cores! please specifiy the num_cores!")
    
    # 1. to find the center sequence
    print("-----------RUN-----------")
    print("cores:", num_cores)
    print("Loading", end="")
    idxC = findCenterSeq(strs, num_cores)
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