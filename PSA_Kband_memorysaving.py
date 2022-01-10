'''
pairwise sequence alignment 
- Affine gap penalty
- Kband
- memory saving

author: Juntao Chen
date: 12.27.2021
'''

match = 1
mismatch = -2
d = 3
e = 1


def score(xi:str , yi:str):
    """return the score of match or mismatch

    Args:
        xi: char one
        yi: char twe
    Returns:
        score
    """
    return match if xi == yi else mismatch


def InsiderStrip(i:int, j:int, k:int, diff=0):
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
    """
    Init the dynamic programming matrix, byte matrix
    """
    bt = [[chr(0)]*(diff+2*k+1) for _ in range(m+1)]
    pm = [[-float('Inf')]*(diff+2*k+1) for _ in range(3)]
    pm[0][k] = 0
    bt[0][k] = chr(16)
    for j in range(1, diff + k + 1):
        pm[1][j+k] = -d - e * (j-1)
        bt[0][j+k] = chr(8)
    for i in range(1, k + 1):
        bt[i][k-i] = chr(3)
    return bt, pm


def InitTwo(i:int, k:int, diff:int):
    """
    Init the dynamic programming matrix
    """
    pm2 = [[-float('Inf')]*(diff+2*k+1) for _ in range(3)]
    if (i < k+1):
        pm2[2][k-i] = -d - e * (i - 1)
    return pm2


def ChooseWay(p0:float, p1:float, p2:float, state=True):
    """
    choose the trace path or state
    
    when state true, return state value

    state false, return path value 
    """
    if p0 >= p1:
        if p0 >= p2:
            return 16 if state else 0
        else:
            return 48 if state else 2
    elif p1 >= p2:
        return 32 if state else 1
    else:
        return 48 if state else 2


def TraceBack(pm, A:str, B:str, k:int, channel:int):
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
    i = m
    b_j = n
    j = n - m + k

    # to get the aligned seqs
    while (i > 0 or j > k):
        if channel == 0:
            channel = parse(pm[i][j], 0)
            seq_A += A[i-1]
            seq_B += B[b_j-1]
            i -= 1
            b_j -= 1

        elif channel == 1:
            channel = parse(pm[i][j], 1)
            seq_A += '-'
            seq_B += B[b_j-1]
            b_j -= 1
            j -= 1

        elif channel == 2:
            channel = parse(pm[i][j], 2)
            seq_A += A[i-1]
            seq_B += '-'
            i -= 1
            j += 1

        else:
            print(i, j, channel)
            raise ValueError("wrong channel = " + channel)

    return seq_A[::-1], seq_B[::-1]


def parse(b:str, s:int):
    """
    compute the path value

    Args:
        b: byte value
        s: path value

    Return:
        next step path value 
    """
    b = ord(b)
    b = (b >> (4 - s * 2))
    return (b & 3) - 1


def PSA_AGP_Kband(A:str, B:str, ms=1, mis=-2, opengap=3, exgap=1, get_score=0):
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
    match = ms
    mismatch = mis
    d = opengap
    e = exgap
    
    # len(A)=0 or len(B)=0
    if len(A) == len(B) == 0:
        return 0, '', ''
    elif len(A) == 0:
        return -e * len(B) - d + e, '-'*len(B), B
    elif len(B) == 0:
        return -e * len(A) - d + e, A, '-'*len(A)
    # n>=m
    # record the loc of A , B
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
        bt, pm = Init(m, k, diff)
        for i in range(1, m+1):
            pm2 = InitTwo(i, k, diff)
            for _ in range(-k, diff+k+1):
                j = _
                if 1 <= j + i <= n:
                    j += k
                    bt1 = bt2 = bt3 = 0
                    # t : A[i] ~ B[j]
                    bt1 = ChooseWay(pm[0][j], pm[1][j], pm[2][j])
                    pm2[0][j] = max(pm[0][j], pm[1][j], pm[2][j]) + score(A[i-1], B[j+i-k-1])

                    if InsiderStrip(i, j+i-k-1, k, diff):
                        # x : B[j] ~ _
                        pm2[1][j] = max(pm2[0][j-1] - d, pm2[1][j-1] - e)
                        bt2 = 4 if pm2[0][j-1] - d > pm2[1][j-1] - e else 8

                    if InsiderStrip(i-1, j+i-k, k, diff):
                        # y : A[i] ~ _
                        pm2[2][j] = max(pm[0][j+1] - d, pm[2][j+1] - e)
                        bt3 = 1 if pm[0][j+1] - d > pm[2][j+1] - e else 3

                    bt[i][j] = chr(bt1 + bt2 + bt3)

            pm = pm2
        new = max(pm[0][diff+k], pm[1][diff+k], pm[2][diff+k])
        if old == new or (k * 2) > m:
            if get_score:
                return new, k
            else:
                break
        else:
            old = new
            k *= 2

    channel = ChooseWay(pm[0][diff+k], pm[1][diff+k], pm[2][diff+k], False)
    SeqA, SeqB = TraceBack(bt, A, B, k, channel)

    # exchange the loc of A & B
    if state_ex:
        SeqA, SeqB = SeqB, SeqA

    return new, SeqA, SeqB


if __name__ == "__main__":
    A = "ACTGACGTAA"
    B = "ACGTGCATTAG"

    print(PSA_AGP_Kband(A, B))
