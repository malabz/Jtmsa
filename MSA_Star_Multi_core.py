'''
star alignment with multiprocessing
author: Juntao Chen
'''

import multiprocessing as mp
import numpy as np
from PSA_Kband import *

# to find the center sequence
def find_censeq(name, param):
    
    S = param[0]
    k = param[1]
    nums = param[2]

    s_psa = [[None]*len(S) for _ in range(len(S))]

    for i in range(len(S)):
        for j in range(len(S)):
            if i < j:
                if (i + j) % nums == k:
                    tmp, _, _ = PSA_AGP_Kband(S[i], S[j])
                    s_psa[i][j] = s_psa[j][i] = tmp

    return s_psa

def insertGap(mark, seq):
    res = ""
    length = len(mark)
    for i in range(length):
        res += "-" * mark[i]
        if i < length - 1: 
            res += seq[i]
    return res

def MSA_star_Multicores(S):

    num_cores = mp.cpu_count()//2
    pool = mp.Pool(processes = num_cores)
    s_psa = [[0]*len(S) for _ in range(len(S))]

    print("-----------RUN-----------")
    print("cores:", num_cores)
    print("Loading", end = "")
    # 1. to find the center sequence
    param_dict = {
        'task1': (S, 0, num_cores),
        'task2': (S, 1, num_cores),
        'task3': (S, 2, num_cores),
        'task4': (S, 3, num_cores),
        'task5': (S, 4, num_cores),
        'task6': (S, 5, num_cores),
        'task7': (S, 6, num_cores),
        'task8': (S, 7, num_cores)
        }
    pool = mp.Pool(num_cores)
    results = [pool.apply_async(find_censeq, args=(name, param)) for name, param in param_dict.items()]
    results = [p.get() for p in results]
    print("=>>>>>")
    for i in range(len(S)):
        for j in range(len(S)):
            if i < j:
                ij = 0
                state = 0
                for k in range(num_cores):
                    if results[k][i][j]:
                        if not state:
                            ij = results[k][i][j]
                        else:
                            raise ValueError("ci wrong")
                s_psa[i][j] = s_psa[j][i] = ij
    C = np.argmax(np.sum(s_psa, axis=0))
    
    print("Loaded")
    print("center seq:", ','.join(S[C]))

    # 2. do pairwise alignments
    Strings = []
    for i in range(len(S)):
        if i != C:
            _, tmp1, tmp2 = PSA_AGP_Kband(S[C], S[i])
            Strings.append([tmp1, tmp2])

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