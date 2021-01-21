'''' 
Levenshtein distance
input: sequence A and B
output: distance d
author: Juntao Chen
'''


def LSD(A, B):

    m = len(A)
    n = len(B)

    p = [[0]*(n+1) for _ in range(m+1)]
    # s = [['']*(n+1) for _ in range(m+1)]

    for i in range(m+1):
        p[i][0] = i
    for j in range(n+1):
        p[0][j] = j

    for i in range(0,m):
        for j in range(0,n):
            if A[i] == B[j]:
                p[i+1][j+1] = p[i][j]
                # s[i+1][j+1] = '='
            else:
                p_min = min(p[i+1][j], p[i][j+1], p[i][j])
                p[i+1][j+1] = p_min + 1

                # if p[i+1][j] == p_min:
                #     s[i+1][j+1] = 'i'
                # if p[i][j+1] == p_min:
                #     if s[i+1][j+1] != '':
                #         s[i+1][j+1] = s[i+1][j+1] + ',d'
                #     else:
                #         s[i+1][j+1] = 'd'
                # if p[i][j] == p_min:
                #     if s[i+1][j+1] != '':
                #         s[i+1][j+1] = s[i+1][j+1] + ',r'
                #     else:
                #         s[i+1][j+1] = 'r'
    return p[-1][-1]


if __name__ == "__main__":
    A = ''
    B = ''
    d = LSD(A, B)
    print(d)



