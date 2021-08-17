'''
pair sequence alignment
Linear gap penalty
Affine gap penalty
Kband
input A and B
output SP_score, align_A and align_B
author: Juntao Chen
'''

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


# Affine gap penalty ~ PSA ~ Kband
def PSA_AGP_Kband(A, B, d=3, e=1, get_score = 0):
    '''len(A)=0 or len(B)=0'''
    if len(A) == len(B) == 0:
        return 0, '', ''
    elif len(A) == 0:
        return -2*len(B), '-'*len(B), B
    elif len(B) == 0:
        return -2*len(A), A, '-'*len(A)
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
                k //= 2
                break
    
    i = m
    b_j = n
    j = diff + k
    seq_A = ""
    seq_B = ""
    score = score_max(t[i][j], x[i][j], y[i][j])
    # print(t)
    # print(x)
    # print(y)

    # to get the aligned seqs
    while (i > 0 or j > k):
        # print(score, i, j, b_j, k)
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
            raise ValueError('Wrong!')
    
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
        if s1[i] != "-" and s2[i] != "-":
            if gap1:
               gap1 = 0
            if s1[i] == s2[i]:
                score_two += m
            else:
                score_two += mis
        elif s1[i] != "-" or s2[i] != "-":
            if gap1 == 0:
                score_two -= d
                gap1 = 1
            else:
                score_two -= e
    return score_two

if __name__ == "__main__":
    A = "AAACCCAAACCCAAEEACCCEE"
    B = "ABCDAAACCCAAAEECCC"
    print(PSA_AGP_Kband(A, B))
    