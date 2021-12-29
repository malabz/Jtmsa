'''
pairwise sequence alignment 
- linear gap penalty
- divide conquer

space complexity: O(n)
time complexity: O(mnlogn)

author: Juntao Chen
date: 12.29.2021
'''

class PSA_DC(object):
    """
    pairwise sequence alignment with divide conquer method

    Args:
        A: sequence 1
        B: sequence 2
        alA: aligned seq 1
        alB: aligned seq 2
    """
    def __init__(self, A, B):
        self.g = -1
        self.A = A
        self.B = B
        # record the match location in t or gap
        self.Align_s = ['']*len(A)
        self.traceback()
        
    # mathch 1 mismatch -1
    def p(self, i, j):
        return int(self.A[i]==self.B[j]) - int(self.A[i]!=self.B[j])

    def match(self, si, tj):
        return int(si==tj) - int(si!=tj)

    def BestScore(self, s1, t1):
        """
        compute the sp score of two seqs
        """
        m = len(s1)
        n = len(t1)
        a = [0] * (n+1)
        
        for j in range(0, n+1):
            a[j] = j * self.g

        for i in range(1, m+1):
            old = a[0]
            a[0] = i * self.g
            for j in range(1, n+1):
                temp = a[j]
                a[j] = max(a[j] + self.g, old + self.match(s1[i-1], t1[j-1]), a[j-1] + self.g)
                old = temp
        return a

    # Divide and conquer strategy
    def Align(self, a, b, c, d):
        if self.A[a:b] == '' or self.B[c:d] == '':
            if self.A[a:b]:
                self.Align_s[a:b] = ['-']*len(self.A[a:b])
            return
        elif '' not in self.Align_s:
            return
        else:
            i = (a+b)//2
            # dimension: d-c
            pref_sim = self.BestScore(self.A[a:i], self.B[c:d])
            t_rev = self.B[c:d][::-1]
            s_rev = self.A[i+1:b][::-1]
            suff_sim = self.BestScore(s_rev, t_rev)
            suff_sim = suff_sim[::-1]

            posmax = c - 1
            typemax = '-'
            Vmax = pref_sim[0] + self.g + suff_sim[0]

            # find the best match (i,j) i from s, j from t
            for j in range(c, d):
                if pref_sim[j-c] + self.p(i,j) + suff_sim[j+1-c] > Vmax:
                    posmax = j
                    typemax = 'N'
                    Vmax = pref_sim[j-c] + self.p(i,j) + suff_sim[j-c+1]

                if pref_sim[j-c+1] + self.g + suff_sim[j-c+1] > Vmax:
                    posmax = j
                    typemax = '-'
                    Vmax = pref_sim[j-c] + self.g + suff_sim[j-c]
            # i match gap
            if typemax == '-':
                self.Align_s[i] = '-' 
                self.Align(a, i, c, posmax+1)
                self.Align(i+1, b, posmax+1, d)
            # i match j
            else:
                self.Align_s[i] = posmax
                self.Align(a, i, c, posmax)
                self.Align(i+1, b, posmax+1, d)

    def traceback(self):
        # compute the record matrix
        self.Align(0, len(self.A), 0, len(self.B))
        idxB_old = -1
        self.alA = ""
        self.alB = ""
        for i in range(len(self.Align_s)):
            if self.Align_s[i] == '-':
                self.alA += self.A[i]
                self.alB += '-'
            else:
                for j in range(idxB_old + 1, self.Align_s[i]):
                    self.alA += '-'
                    self.alB += self.B[j]
                self.alA += self.A[i]
                self.alB += self.B[self.Align_s[i]]
                idxB_old = self.Align_s[i]
        for i in range(idxB_old + 1, len(self.B)):
            self.alA += '-'
            self.alB += self.B[i]
            

if __name__ == "__main__":
    # input sequences
    B = 'ACG'
    A = 'TTT'

    psa = PSA_DC(A, B)

    print(psa.alA)
    print(psa.alB)