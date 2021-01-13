# pair sequence alignment
# Linear gap penalty
# Affine gap penalty
# Kband
# author: Juntao Chen

import os
import datetime
from pprint import pprint
from Extract_data import read_fasta

# match 1
# mismatch -2
# gap
# d = 3
# e = 1

# return the score of match or mismatch
def s(xi,yi,m=1,mis=-2):
    if xi == yi:
        return m
    else:
        return mis


# to record the direction
def score_max(t,x,y):
    m = max(t, x, y)
    if t == m:
        return 't'
    elif x == m:
        return 'x'
    elif y == m:
        return 'y'


# to return the existence of the index ~ for kband function
def InsiderStrip(i, j, k, diff=0):
    return (-k <= j - i <= k + diff)


# PSA ~ linear gap penalty
def PSA_LGP(A, B, g = -1):
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A)>len(B):
        A, B = B, A
        state_ex = 1
    n = len(B)
    m = len(A)

    p = [[0]*(n+1) for _ in range(m+1)]

    # init
    for i in range(1, m+1):
        p[i][0] = g*i
    for j in range(1, n+1):
        p[0][j] = g*j

    for i in range(0, m):
        for j in range(0, n):
            p[i+1][j+1] = max(p[i][j] + s(A[i], B[j]), p[i][j+1] + g, p[i+1][j] + g)

    i = m
    j = n
    seq_A = ""
    seq_B = ""
    
    while (i > 0 or j > 0):
        if i > 0 and j > 0:
            if p[i][j] == p[i-1][j-1] + s(A[i-1], B[j-1]):
                seq_A += A[i-1]
                seq_B += B[j-1]        
                i -= 1
                j -= 1
                continue
        if i > 0:
            if p[i][j] == p[i-1][j] + g:
                seq_A += A[i-1]
                seq_B += '-'        
                i -= 1
                continue
        if j > 0:
            if p[i][j] == p[i][j-1] + g:
                seq_A += '-'
                seq_B += B[j-1]        
                j -= 1
                continue
        else:
            print(i,j)
            raise RuntimeError('Error')

    # exchange the loc of A & B
    if state_ex:
        seq_A, seq_B = seq_B, seq_A
    
    return p[-1][-1], seq_A[::-1], seq_B[::-1]


# PSA ~ Affine gap penalty
def PSA_AGP(A, B, d=3, e=1):
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A)>len(B):
        A, B = B, A
        state_ex = 1
    n = len(B)
    m = len(A)

    t = [[-float('Inf')]*(n+1) for _ in range(m+1)]
    x = [[-float('Inf')]*(n+1) for _ in range(m+1)]
    y = [[-float('Inf')]*(n+1) for _ in range(m+1)]

    # init
    t[0][0] = 0
    for i in range(1, m+1):
        y[i][0] = -d - e*(i-1)

    for j in range(1, n+1):
        x[0][j] = -d - e*(j-1)

    for i in range(1, m+1):
        for j in range(1, n+1):
            # x :  _ ~ B[j]
            x[i][j] = max(t[i][j-1]-d, x[i][j-1]-e)
            # y : A[i] ~ _ 
            y[i][j] = max(t[i-1][j]-d, y[i-1][j]-e)
            # t : A[i] ~ B[j]
            t[i][j] = max(t[i-1][j-1], x[i-1][j-1], y[i-1][j-1]) + s(A[i-1], B[j-1])

    i = m
    j = n
    seq_A = ""
    seq_B = ""

    score_ = max(t[-1][-1], x[-1][-1], y[-1][-1])
    score = score_max(t[i][j], x[i][j], y[i][j])

    while (i > 0 or j > 0):
        if score == 't' and i>0 and j>0:
            if t[i][j] == t[i-1][j-1] + s(A[i-1], B[j-1]) and i>1 and j>1:
                score = 't'
            elif t[i][j] == x[i-1][j-1] + s(A[i-1], B[j-1]) and j>1:
                score = 'x'
            elif t[i][j] == y[i-1][j-1] + s(A[i-1], B[j-1]) and i>1:
                score = 'y'
            seq_A += A[i-1]
            seq_B += B[j-1]
            i -= 1
            j -= 1
        elif score == 'x' and j>0:
            if x[i][j] == t[i][j-1] - d and i>=1 and j>1:
                score = 't'
            elif x[i][j] == x[i][j-1] - e:
                score = 'x'
            seq_A += '-'
            seq_B += B[j-1]
            j -= 1
        elif score == 'y' and i>0:
            if y[i][j] == t[i-1][j] - d and i>1 and j>=1:
                score = 't'
            elif y[i][j] == y[i-1][j] - e:
                score = 'y'
            seq_A += A[i-1]
            seq_B += '-'
            i -= 1

    # exchange the loc of A & B
    if state_ex:
        seq_A, seq_B = seq_B, seq_A

    return score_, seq_A[::-1], seq_B[::-1]


# PSA ~ Kband ~ linear penalty
def PSA_LGP_Kband(A, B, g=-1):
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A)>len(B):
        A, B = B, A
        state_ex = 1
    n = len(B)
    m = len(A)
    diff = n - m
    k = 1
    p = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    old = -float('Inf')

    while k <= m:
        # init
        for i in range(k+1):
            p[i][k-i] = i * g
        for j in range(1, k+1+diff):
            p[0][j+k] = j * g
        p[0][k] = 0
        
        for i in range(1, m+1):
            for d in range(-k, k+1+diff):
                j = d
                if 1 <= j + i <= n:
                    j += k
                    p[i][j] = p[i-1][j] + s(A[i-1],B[j+i-k-1])
                    if InsiderStrip(i-1, j+i-k, k, diff):
                        p[i][j] = max(p[i][j], p[i-1][j+1] + g)
                    if InsiderStrip(i, j-1+i-k, k, diff):
                        p[i][j] = max(p[i][j], p[i][j-1] + g)
        if old == p[-1][-1-k]:
            # pprint(a)
            break
        else:
            old = p[-1][-1-k]
            k *= 2
            if k <= m:
                p = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
            else:
                k /= 2
                break
    i = m
    b_j = n
    j = diff + k

    seq_A = ""
    seq_B = ""
    
    # to get the aligned seqs
    while (i > 0 or j > k):
        if i > 0 and j >= 0:
            if p[i][j] == p[i-1][j] + s(A[i-1], B[b_j-1]):
                seq_A += A[i-1]
                seq_B += B[b_j-1]        
                i -= 1
                b_j -= 1
                continue

        if i > 0 and j + 1 <= 2 * k + diff:
            if p[i][j] == p[i-1][j+1] + g:
                seq_A += A[i-1]
                seq_B += '-'        
                i -= 1
                j += 1
                continue
        
        if j > 0:
            if p[i][j] == p[i][j-1] + g:
                seq_A += '-'
                seq_B += B[b_j-1]        
                b_j -= 1
                j -= 1
                continue
        else:
            print(i,j,b_j)
            raise RuntimeError('Error')


    # exchange the loc of A & B
    if state_ex:
        seq_A, seq_B = seq_B, seq_A
        
    return p[-1][-1-k], seq_A[::-1], seq_B[::-1]


# Affine gap penalty ~ PSA ~ Kband
def PSA_AGP_Kband(A, B, d=3, e=1, get_score = 0):
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A)>len(B):
        A, B = B, A
        state_ex = 1
    n = len(B)
    m = len(A)
    diff = n - m
    k = 1

    t = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    x = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    y = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    old = -float('Inf')

    # to compute the optimal score
    while k <= m:
        # init
        t[0][k] = 0
        for i in range(1, k+1):
            y[i][k-i] = -d - e*(i-1)
        for j in range(1, k+1+diff):
            x[0][j+k] = -d - e*(j-1)

        for i in range(1, m+1):
            for _ in range(-k, diff+k+1):
                j = _
                if 1 <= j + i <= n:
                    j += k
                    # t : A[i] ~ B[j]
                    t[i][j] = max(t[i-1][j], x[i-1][j], y[i-1][j]) + s(A[i-1], B[j+i-k-1])
                    
                    if InsiderStrip(i, j+i-k-1, k, diff):
                        # x : B[j] ~ _ 
                        x[i][j] = max(t[i][j-1]-d, x[i][j-1]-e)
                    
                    if InsiderStrip(i-1, j+i-k, k, diff):
                        # y : A[i] ~ _ 
                        y[i][j] = max(t[i-1][j+1]-d, y[i-1][j+1]-e)
        
        if old == max(t[-1][-1-k], x[-1][-1-k], y[-1][-1-k]):
            if get_score:
                return old,k
            else:
                break
        else:
            old = max(t[-1][-1-k], x[-1][-1-k], y[-1][-1-k])
            k *= 2
            if k <= m:
                    t = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
                    x = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
                    y = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
            else:
                k /= 2
                break
    
    i = m
    b_j = n
    j = diff + k
    seq_A = ""
    seq_B = ""
    score = score_max(t[i][j-k], x[i][j-k], y[i][j-k])

    # to get the aligned seqs
    while (i > 0 or j > k):
        if score == 't' and i > 0 and j >= 0:
            if t[i][j] == t[i-1][j] + s(A[i-1], B[b_j-1]) and i > 1 and j >= 0:
                score = 't'
            elif t[i][j] == y[i-1][j] + s(A[i-1], B[b_j-1]) and i > 1:
                score = 'y'
            elif t[i][j] == x[i-1][j] + s(A[i-1], B[b_j-1]) and j > 0:
                score = 'x'
            seq_A += A[i-1]
            seq_B += B[b_j-1]
            i -= 1
            b_j -= 1

        elif score == 'x' and j > 0:
            if x[i][j] == x[i][j-1] - e:
                score = 'x'
            elif x[i][j] == t[i][j-1] - d and i >= 1:
                score = 't'

            seq_A += '-'
            seq_B += B[b_j-1]
            b_j -= 1
            j -= 1

        elif score == 'y' and i > 0 and j + 1 <= 2 * k + diff:  
            if y[i][j] == y[i-1][j+1] - e:
                score = 'y'
            elif y[i][j] == t[i-1][j+1] - d and i > 1 and j >= 0:
                score = 't'

            seq_A += A[i-1]
            seq_B += '-'
            i -= 1
            j += 1

        else:
            print(i,j,score)
            raise RuntimeError('Wrong!')
    
    # exchange the loc of A & B
    if state_ex:
        seq_A, seq_B = seq_B, seq_A
    
    return max(t[-1][-1-k], x[-1][-1-k], y[-1][-1-k]), seq_A[::-1], seq_B[::-1]


# compute the SP of two seqs
def Compute_two(s1, s2, d=3, e=1, m=1, mis=-2):
    if len(s1) != len(s2):
        print(s1, s2)
        raise ValueError("the length of s1 and s2 is wrong!")
    score_two = 0
    gap1 = 0
    for i in range(len(s1)):
        if s1[i] != "_" and s2[i] != "_":
            if gap1:
               gap1 = 0
            if s1[i] == s2[i]:
                score_two += m
            else:
                score_two += mis
        elif s1[i] != "_" or s2[i] != "_":
            if gap1 == 0:
                score_two -= d
                gap1 = 1
            else:
                score_two -= e
    return score_two


if __name__ == "__main__":

    name = ['SARS','dog_eye','Homo_sapiens','SCML4','genome']
    for i in range(len(name)):
        # i = 3
        begin_time = datetime.datetime.now()
        
        strs = read_fasta(filename=name[i]+'_aligned.fasta')
        # strs = read_fasta(filename=name[i]+'.fasta')[0:2]
        file_path = os.getcwd() + "/data/" + name[i] + "_aligned.fasta"

        print(len(strs))
        for j in range(1,2):
            
            # A = strs[0]
            # B = strs[j]
            # print(name[i],":",len(A),len(B))

            # _, A_aligned, B_aligned = PSA_AGP_Kband(A, B, d = 3, e = 1)

            end_time = datetime.datetime.now()
            
            score_1 = Compute_two(strs[0], strs[1])
            score_2 = Compute_two(strs[2], strs[3])
            # score_3 = Compute_two(strs[4], strs[5])
            # score_4 = Compute_two(strs[6], strs[7])

            with open(file_path, 'a+') as f:
                f.write('\n' + 'SP score of Kband:' + str(score_1))
                f.write('\n' + 'SP score of MAFFT:' + str(score_2))
                # f.write('\n' + 'SP score of Kband:' + str(score_3))
                # f.write('\n' + 'SP score of MAFFT:' + str(score_4))

            # with open(file_path, 'w') as f:
            #     f.write('\n' + '> 0 ' + 'Kband  Time cost:' + str(end_time-begin_time) + ' length: '+ str(len(A)) + '\n')
            #     f.write(A_aligned)
            #     f.write('\n' + '> ' + str(j) + ' length: '+ str(len(B)) + '\n')
            #     f.write(B_aligned)
            #     f.write('\n')

        print((end_time-begin_time).seconds)