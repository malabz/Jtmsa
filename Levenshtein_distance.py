'''' 
Levenshtein distance
space complexity: O(n)
time complexity: O(mn)

author: Juntao Chen
date: 12.28.2021
'''


def LSD(A, B):
    """
    compute Levenshtein distance

    Args:
        A: string / list 
        B: string / list

    Returns:
        distance
    """
    m = len(A)
    n = len(B)
    p = [0] * (n + 1)
    for j in range(n+1):
        p[j] = j

    for i in range(1, m + 1):
        old = p[0]
        p[0] = i
        for j in range(1, n + 1):
            temp = p[j]
            if A[i-1] == B[j-1]:
                p[j] = old
            else:
                p[j] = min(p[j-1], old, p[j]) + 1
            old = temp
    return p[-1]


if __name__ == "__main__":
    A = 'AAAA'
    B = 'ABACCAA'
    d = LSD(A, B)
    print(d)