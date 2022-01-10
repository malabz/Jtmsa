'''
pairwise sequence alignment
- Affine gap penalty
- Kband

author: Juntao Chen
date: 12.27.2021
'''

match = 1
mismatch = -2
d = 3
e = 1


def score(xi: str, yi: str):
    """return the score of match or mismatch

    Args:
        xi: char one
        yi: char twe
    Returns:
        score
    """
    return match if xi == yi else mismatch


def ChooseWay(p0: float, p1: float, p2: float):
    """
    choose the trace path
    """
    if p0 >= p1:
        if p0 >= p2:
            return 't'
        else:
            return 'y'
    elif p1 >= p2:
        return 'x'
    else:
        return 'y'


def InsiderStrip(i: int, j: int, k: int, diff=0):
    """judge whether the loc is in the k-band area

    Args:
        i: loc in y
        j: loc in x
        k: kband size
        diff: Sequence length difference
    Returns:
        True or False
    """
    return (-k <= j - i <= k + diff)


def Init(m:int, k:int, diff:int):
    t = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    x = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    y = [[-float('Inf')]*(diff+2*k+1) for _ in range(m+1)]
    # init
    t[0][k] = 0
    for i in range(1, k+1):
        y[i][k-i] = -d - e*(i-1)
    for j in range(1, k+1+diff):
        x[0][j+k] = -d - e*(j-1)
    
    return t, x, y


def TraceBack(pm:list, A:str, B:str, k:int, channel:int):
    """
    Trace back the optimal path, and get the optimal alignment

    Args:
        pm: path matrix
        A: seq 1
        B: seq 2
        k: kband size
        channel: last step's path

    Returns:
        aligned A, aligned B 
    """
    seq_A = ""
    seq_B = ""
    m = len(A)
    n = len(B)
    diff = n - m
    i = m
    b_j = n
    j = diff + k
    t = pm[0]
    x = pm[1]
    y = pm[2]

    # to get the aligned seqs
    while (i > 0 or j > k):
        if channel == 't' and i > 0 and j >= 0:
            if t[i][j] == t[i-1][j] + score(A[i-1], B[b_j-1]) and i > 1 and j >= 0:
                channel = 't'
            elif t[i][j] == y[i-1][j] + score(A[i-1], B[b_j-1]) and i > 1:
                channel = 'y'
            elif t[i][j] == x[i-1][j] + score(A[i-1], B[b_j-1]) and j > 0:
                channel = 'x'
            seq_A += A[i-1]
            seq_B += B[b_j-1]
            i -= 1
            b_j -= 1

        elif channel == 'x' and j > 0:
            if x[i][j] == x[i][j-1] - e:
                channel = 'x'
            elif x[i][j] == t[i][j-1] - d and i >= 1:
                channel = 't'
            seq_A += '-'
            seq_B += B[b_j-1]
            b_j -= 1
            j -= 1

        elif channel == 'y' and i > 0 and j + 1 <= 2 * k + diff:
            if y[i][j] == y[i-1][j+1] - e:
                channel = 'y'
            elif y[i][j] == t[i-1][j+1] - d and i > 1 and j >= 0:
                channel = 't'
            seq_A += A[i-1]
            seq_B += '-'
            i -= 1
            j += 1

        else:
            print(i, j, channel)
            raise ValueError("wrong channel = " + channel)

    return seq_A[::-1], seq_B[::-1]


def PSA_AGP_Kband(A: str, B: str, m=1, mis=-2, opengap=3, exgap=1, get_score=0):
    """
    Affine gap penalty ~ PSA ~ Kband

    Args:
        A: sequence 1
        B: sequence 2

    Returns:
        value1: align score

        value2: aligned seq A

        value3: aligned seq B
    """
    global match, mismatch, d, e
    match = m
    mismatch = mis
    d = opengap
    e = exgap

    # len(A)=0 or len(B)=0
    if len(A) == len(B) == 0:
        return 0, '', ''
    elif len(A) == 0:
        return -2*len(B), '-'*len(B), B
    elif len(B) == 0:
        return -2*len(A), A, '-'*len(A)
    # n>=m
    # record the loc of A & B
    state_ex = 0
    if len(A) > len(B):
        A, B = B, A
        state_ex = 1

    n = len(B)
    m = len(A)
    diff = n - m
    k = 1

    old = -float('Inf')

    # to compute the optimal score
    while k <= m:
        t,x,y = Init(m, k, diff)
        for i in range(1, m+1):
            for _ in range(-k, diff+k+1):
                j = _
                if 1 <= j + i <= n:
                    j += k
                    # t : A[i] ~ B[j]
                    t[i][j] = max(t[i-1][j], x[i-1][j], y[i-1][j]) 
                    t[i][j] += score(A[i-1], B[j+i-k-1])
                    
                    if InsiderStrip(i, j+i-k-1, k, diff):
                        # x : B[j] ~ _
                        x[i][j] = max(t[i][j-1]-d, x[i][j-1]-e)

                    if InsiderStrip(i-1, j+i-k, k, diff):
                        # y : A[i] ~ _
                        y[i][j] = max(t[i-1][j+1]-d, y[i-1][j+1]-e)

        new = max(t[-1][-1-k], x[-1][-1-k], y[-1][-1-k])
        if old == new or (k * 2) > m:
            if get_score:
                return new, k
            else:
                break
        else:
            old = new
            k *= 2

    channel = ChooseWay(t[-1][-1-k], x[-1][-1-k], y[-1][-1-k])
    SeqA, SeqB = TraceBack([t,x,y], A, B, k, channel)

    # exchange the loc of A & B
    if state_ex:
        SeqA, SeqB = SeqB, SeqA

    return new, SeqA, SeqB


if __name__ == "__main__":
    A = "ACTGACGTAA"
    B = "ACGTGCATTAG"

    print(PSA_AGP_Kband(A, B))