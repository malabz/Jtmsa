'''
pairwise sequence alignment with suffix tree

author: Juntao Chen
date: 12.29.2021
'''

import numpy as np

from PSA_Kband import PSA_AGP_Kband
from suffixTree import suffixTree, Node
from score import spscore

class PSA_STree(object):
    def __init__(self, strings:str):
        self.T = strings + '#'
        self.len = len(strings)
        self.root = suffixTree(strings).root
        # set the minimum length of common strings
        self.thorehold = len(strings) // 100 if len(strings) // 100 >= 15 else 15;

    def _edge_length(self, node:Node):
        """
        return edgs's length
        """
        return node.end - node.start + 1
    
    def __walk_down_fcs(self, node:Node, step:int):
        """
        return whether the step excess this edge
        """
        return self._edge_length(node) <= step

    def _dfs_leaves(self, node:Node, results:list, length:int):
        """
        use DPS (depth first search) to search the leaves and return all the index

        Args:
            node: Node
            results: use to save the idx
            length: walk step
        """
        sons = node.children.values()
        for son in sons:
            if son.leaf:
                results.append(son.start-length-self._edge_length(node))
            else:
                tmp = length + self._edge_length(son)
                self._dfs_leaves(son, results, tmp)
        return results

    def _select_prefix(self, s:str):
        """
        to find the longest common strings with the prefix of s

        Returns:
            starts: lcs's index (list)
            length: lcs's length
        """
        node = self.root
        starts = []
        length, step, tag = 0, 0, 0
        while node.children.get(s[step+length]) is not None:
            node = node.children[s[step+length]]
            step += 1
            if len(s) >= step + length:
                if self.__walk_down_fcs(node, step):
                    length += step
                    step = 0
                    if len(s) == length:
                        break
                    else:
                        continue
                if len(s) == step + length:
                    break
            else:
                break

            if self.T[node.start + step] != s[step+length]:
                break

            while self.T[node.start + step] == s[step+length]:
                step += 1
                if len(s) >= step + length:
                    if self.__walk_down_fcs(node, step):
                        length += step
                        step = 0
                        if len(s) == length:
                            tag = 1
                        break
                    if len(s) == step + length:
                        tag = 1
                        break
                else:
                    tag = 1
                    break
            if tag:
                break
            if self.T[node.start + step] != s[step+length]:
                break
        if node.leaf:
            starts.append(node.start-length)
        elif length + step >= self.thorehold:
            starts = self._dfs_leaves(node, starts, length)
        else:
            starts = []
        return starts, length + step

    def _findCommonStrings(self, s2:str):
        """
        find the all common strings between s2 and T

        Returns:
            results: [[idxB, length, idxA's list]]
        """
        # results;
        index = 0
        results = []
        while index <= len(s2) - 1:
            starts, length = self._select_prefix(s2[index:])
            if starts and length > 1:
                results.append([index, length, starts])
                index += length
            else:
                index += 1
        del_k = [k for k, i in enumerate(results) if i[2] == []]
        return [results[i] for i in range(len(results)) if i not in del_k]
    
    def _multi_To_one(self, results:list):
        """
        choose a closest idxA for each idxB and delete the length smaller than thorehold

        Args:
            results: return by function _findCommonStrings()
        """
        start = float('Inf')
        end = -float('Inf')
        length = results[-1][0] + results[-1][1]
        for result in results:
            s = min(result[2])
            e = max(result[2]) + result[1]
            start = s if s < start else start
            end = e if e > end else end
        delete_k = []
        for k, result in enumerate(results):
            if result[1] < self.thorehold: 
                delete_k.append(k)
            elif len(result[2]) > 1:
                temp = np.array(result[2])
                temp = (temp / (end - start)) - result[0] / length
                temp = list(abs(temp))
                index = temp.index(min(temp))
                result[2] = result[2][index]
            else:
                result[2] = result[2][0]
        return [results[i] for i in range(len(results)) if i not in delete_k]

    def _trace_back(self, results:list, p:list):
        """
        trace back the index used by function _select_dp()
        """
        i = p.index(max(p))
        track = [i]
        while i > 0:
            j = i - 1
            if p[i] == results[i][1]: break
            while j >= 0:
                if results[i][2] >= (results[j][2]+results[j][1]) and p[i] == p[j] + results[i][1]:
                    track.append(j)
                    i = j
                    break
                j -= 1
        return track[::-1]
    
    def _select_dp(self, results:list):
        """
        select the longest and compatible common strings using DP

        Args:
            results: return by function _multi_To_one()
        """
        m = len(results)
        if m == 0: 
            return [], 0
        if m == 1: 
            return results, results[0][1]/self.len
        p = [i[1] for i in results]
        for i in range(1, m):
            for j in range(i):
                if results[i][2] >= (results[j][2]+results[j][1]):
                    p[i] = max(p[i], p[j] + results[i][1])
        return [results[i] for i in self._trace_back(results, p)], max(p)/self.len

    def _select_CommonStrings(self, s2:str):
        """
        select the compatible common strings between s2 and T
        """
        results = self._findCommonStrings(s2)
        if results == []: return [], 0
        return self._select_dp(self._multi_To_one(results))

    def align(self, s2:str):
        results, sim = self._select_CommonStrings(s2)
        A, B = [], []
        align_A, align_B = [], []
        s_a, s_b = 0, 0
        for r in results:
            A.append(self.T[s_a:r[2]])
            B.append(s2[s_b:r[0]])
            s_a = r[1] + r[2]
            s_b = r[1] + r[0]
        A.append(self.T[s_a:-1])
        B.append(s2[s_b:])
        for i in range(len(A)):
            _, temp_A, temp_B = PSA_AGP_Kband(A[i], B[i])
            align_A.append(temp_A)
            align_B.append(temp_B)
        s_A, s_B = '', ''
        for j in range(len(results)):
            s_A += align_A[j]
            s_B += align_B[j]
            s_A += s2[results[j][0]:(results[j][0]+results[j][1])]
            s_B += s2[results[j][0]:(results[j][0]+results[j][1])]
        s_A += align_A[-1]
        s_B += align_B[-1]
        score = spscore([s_A, s_B])
        return sim, score, s_A, s_B

if __name__ == "__main__":
    A = "AACTCTAAACTCTAAACTCTAAACTCTAAACTCTAAACTCTAAACTCTA"
    B = "AACTAAACTAAACTAAACTAAACTAAAAAACTCTAAACTCTAAACTCTA"
    stree = PSA_STree(A)
    print(stree.align(B))

