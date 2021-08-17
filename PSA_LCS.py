'''
PSA with LCS longest common sequence (Nakatsu LCS)
self.align()
input A and B
output spscore, align_A, align_B
author: Juntao Chen
'''

from PSA_Kband import Compute_two, PSA_AGP_Kband

class PSA_LCS(object):
    def __init__(self, A:str, B:str):
        self.state_loc = 0
        if len(A) > len(B):
            A, B = B, A
            self.state_loc = 1
        self.A = A
        self.B = B
        self.m = len(A)
        self.n = len(B)
        self.length = -1
    
    def _min_j(self, l1, l2, i):
        if l1 == -1:
            return -1
        if l2 == -1:
            l2 = self.n + 1
        for j in range(l1+1, l2):
            if self.A[i-1] == self.B[j-1]:
                return j
        return -1
    
    def _trace_back(self, L, k, i):
        track = []
        while k > 0 and i > 0:
            if L[k][i-1] == -1:
                L[k][i-1] = float('Inf')
            if L[k-1][i-1] < L[k][i] < L[k][i-1]:
                track.append([k, L[k][i], i])
                k -= 1
                i -= 1
            else:
                i -= 1
        track = track[::-1]
        # print(track, len(track))
        return self._merge(track)

    def _judge(self, t1, t2):
        return True if t1[1] + 1 == t2[1] and t1[2] + 1 ==t2[2] else False

    def _merge(self, track):
        merge_track = []
        start, end = [], []
        i = 0
        for i in range(len(track)-1):
            if self._judge(track[i], track[i+1]):
                if start == []: start = track[i]
                end = track[i+1]
            else:
                if start:
                    temp = [[start[1]-1, end[1]-1], [start[2]-1, end[2]-1]]
                    merge_track.append(temp)
                    start, end = [], []
                else:
                    temp = [[track[i][1]-1,track[i][1]-1], [track[i][2]-1,track[i][2]-1]]
                    merge_track.append(temp)
        if start:
            temp = [[start[1]-1, end[1]-1], [start[2]-1, end[2]-1]]
            merge_track.append(temp)
            start, end = [], []
        else:
            temp = [[track[i][1]-1,track[i][1]-1], [track[i][2]-1,track[i][2]-1]]
            merge_track.append(temp)
        del_len = [k for k, i in enumerate(merge_track) if i[0][1]-i[0][0] <= self.n//100]
        new_track = [merge_track[i] for i in range(len(merge_track)) if i not in del_len]
        return new_track

    def LCS(self):
        L = [[-1]*(self.m+1) for _ in range(self.m+1)]
        for i in range(self.m+1):
            L[0][i] = 0
        i, j, k = 0, 0 ,0
        while j <= self.m:
            i = 1 + j
            k = 1
            while j + k <= self.m:
                temp = self._min_j(L[k-1][i-1], L[k][i-1], i)
                L[k][i] = temp if temp > -1 else L[k][i-1]
                if L[k][i] == -1:
                    break
                i += 1
                k += 1
            if j + k> self.m:
                if L[k-1][i-1] != -1:
                    self.length = k - 1
                    break
            j += 1
        track = self._trace_back(L, k-1, i-1)
        return track
    
    def sim(self, track):
        length = 0
        for i in track: length += (i[0][1]-i[0][0]+1)
        return length/self.n

    def align(self):
        track = self.LCS()
        sim = self.sim(track)
        strs_A, strs_B = [], []
        align_A, align_B = [], []
        s_a, s_b = 0, 0
        for t in track:
            strs_B.append(self.B[s_b:t[0][0]])
            strs_A.append(self.A[s_a:t[1][0]])
            s_a = t[1][1] + 1
            s_b = t[0][1] + 1
        strs_B.append(self.B[s_b:])
        strs_A.append(self.A[s_a:])

        for i in range(len(strs_A)):
            _, temp_A, temp_B = PSA_AGP_Kband(strs_A[i], strs_B[i])
            align_A.append(temp_A)
            align_B.append(temp_B)
        
        s_A, s_B = '', ''
        for j in range(len(track)):
            s_A += align_A[j]
            s_B += align_B[j]
            s_A += self.B[track[j][0][0]:track[j][0][1]+1]
            s_B += self.B[track[j][0][0]:track[j][0][1]+1]
        s_A += align_A[-1]
        s_B += align_B[-1]

        score = Compute_two(s_A, s_B)
        if self.state_loc: s_A, s_B = s_B, s_A
        return sim, score, s_A, s_B

if __name__ == "__main__":
    pass

