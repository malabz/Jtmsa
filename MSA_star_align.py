'''
MSA---star alignment
author: Juntao Chen
'''

from PSA import Compute_two, PSA_AGP_Kband
from multiprocessing import Process, Pool
from Extract_data import extract_data
from pprint import pprint
import numpy as np
import datetime

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

    s = 50
    t = 500
    for i in range(0,1):
        s -= 50
        print(s)
        begin_time = datetime.datetime.now()
        strs = extract_data()
        Max, strs_aligned = MSA_star(strs[s:t])
        end_time = datetime.datetime.now()
        run_time = end_time - begin_time
        print("the cost of time:", run_time)