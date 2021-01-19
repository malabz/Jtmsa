# star alignment with multiprocessing
# author: Juntao Chen

from Extract_data import extract_data, read_fasta
from PSA import PSA_AGP_Kband, Compute_two
from pprint import pprint
import multiprocessing as mp
import numpy as np
import datetime
import os

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
        # 'task5': (S, 4, num_cores),
        # 'task6': (S, 5, num_cores),
        # 'task7': (S, 6, num_cores),
        # 'task8': (S, 7, num_cores)
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
                            raise("ci wrong")
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
    for str in Strings:
        if S_aligned:
            j = 0
            k = 0
            marks = [[],[]]
            while j<len(S_aligned[0]) and k<len(str[0]):
                if S_aligned[0][j] == str[0][k]:
                    j += 1
                    k += 1                 
                elif S_aligned[0][j] == "-":
                    marks[1].append(k)
                    j += 1
                elif str[0][k] == "-":
                    marks[0].append(j)
                    k += 1
                else:
                    j += 1
                    k += 1

            tmp_j = len(S_aligned[0]) - j
            tmp_k = len(str[0]) - k
            
            for mark in marks[0]:
                for i in range(len(S_aligned)):
                    S_aligned[i] = S_aligned[i][0:mark] + "-" + S_aligned[i][mark:]
            for mark in marks[1]:
                str[1] = str[1][0:mark] + "-" + str[1][mark:]
            
            if tmp_j:
                str[1] += tmp_j * "-"
            elif tmp_k:
                for i in range(len(S_aligned)):
                    S_aligned[i] += tmp_k * "-"               
            
            if len(S_aligned[0]) != len(str[1]):
                print(tmp_j, tmp_k)
                print(S_aligned[0])
                print(str[0])
                print(str[1])
                raise("the length of seqs have a problem")            
            
            S_aligned.append(str[1])
        else:
            S_aligned.append(str[0])
            S_aligned.append(str[1])
    
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
    
    s = 0
    e = 1000
    begin_time = datetime.datetime.now()
    # strs = extract_data()
    strs = read_fasta()
    print(len(strs))
    Max, strs_aligned = MSA_star_Multicores(strs[s:e])
    end_time = datetime.datetime.now()
    run_time = end_time - begin_time
    print("the cost of time:", run_time)