# -*- coding: utf-8 -*-
'''
BWT-PSA
时间复杂度 O(nlgn)
空间复杂度 O(n)

author: ZhouTong
time :2021.10.14
'''
import numpy as np

from collections import Counter
from PSA_Kband_memorysaving import PSA_AGP_Kband
from score import spscore
from FASTA import readfasta   #提取或写入DNA

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

class PSA_BWT(object):
    def __init__(self, strings):
        self.T_ = strings[::-1]+"#"  # “#”为字符串结束标志
        self.T = strings+"#"
        self.len = len(strings)
        self.thorehold = len(strings) // 10000 if len(strings) // 10000 >= 15 else 15
        self.BWT()

    def bz_sort(self):#倍增排序求 后缀数组，L
        n = self.len + 1
        x = [None] * n  # int，转为int的原串
        sa = [None] * n  #
        y = [None] * n
        c = [None] * n
        for i in range(n):  # 将字符按优先级转为int
            if self.T_[i] == 'A':
                x[i] = 1
            elif self.T_[i] == 'C':
                x[i] = 2
            elif self.T_[i] == 'G':
                x[i] = 3
            elif self.T_[i] == 'T':
                x[i] = 4
            else:
                x[i] = 0
        for i in range(5):
            c[i] = 0
        for i in range(n):
            c[x[i]] += 1  # ACGT各有多少
        for i in range(1, 5):
            c[i] += c[i - 1]  # F列的索引分界点
        for i in range(n - 1, -1, -1):  # 倒着来，是为了保证在当字符串中有相等的字符串时，默认靠前的字符串更小一些。
            c[x[i]] -= 1
            sa[c[x[i]]] = i  # sa[i] 表示的是按第一关键字排序后的第 i 个后缀在原字符串的起始位置。初始化，未排

        # 下面这层循环中p代表rank值不用的字符串的数量，如果p达到n，那么各个字符串的大小关系就已经明了了。
        # k代表当前待合并的字符串的长度，每次将两个长度为j的字符串合并成一个长度为2 * j的字符串，当然如果包含字符串末尾具体则数值应另当别论，但思想是一样的。
        # m同样代表基数排序的元素的取值范围
        k = 1
        m = 5
        p = 1
        while (p < n):
            # 这一部分是按第二关键字排序
            p = 0
            for i in range(n - k, n):  # 将没有第二关键字的后缀的第二关键字设置为无限小，所以要放在 y[] 数组的最前面。
                y[p] = i
                p += 1
            # 看图，下面一行的第二关键字不为0的部分都是根据上面一行的排序结果得到的，且上一行中只有sa[i] >= j的第sa[i]个字符串
            # （这里以及后面指的“第?个字符串”不是按字典序排名来的，是按照首字符在字符串中的位置来的）
            # 的rank才会作为下一行的第sa[i] - j个字符串的第二关键字，而且显然按sa[i]
            # 的顺序rank[sa[i]]是递增的，因此完成了对剩余的元素的第二关键字的排序。
            for i in range(n):  # sa[i] - k 就相当于按第一关键字排序后的第 i 个后缀在原字符串的下标的第前 k 个位置
                if (sa[i] >= k):
                    y[p] = sa[i] - k
                    p += 1
            # 上面第二关键字排序完成后，y[]存放的是按第二关键字排序的字符串下标

            # 按第一关键字进行的基数排序
            for i in range(m):
                c[i] = 0
            for i in range(n):
                c[x[y[i]]] += 1
            #print(c)  C记录了每个排名有多少个重复的，如果每个排名都不重复，则全1，则排序完成
            for i in range(1, m):
                c[i] += c[i - 1]
            #c每步累加，得到每个排名的分界点，下一步再每步从尾部遍历，索引从分界点递减
            for i in range(n - 1, -1, -1):
                c[x[y[i]]] -= 1
                sa[c[x[y[i]]]] = y[i]

            # 结合两个关键字的分别排序，构建出当前的双关键字排序后的结果，也是下一次长度为2k时按第一关键字排序后的结果。
            # y[]本身具有顺序性，即第二关键字小的在前面，而本身这个排序则是将第一关键字小的排在前面，而不改变第一关键字相同的数之间的相对位置，
            # 综合起来，得到的答案便是按双关键字排序后的结果。

            # 下面两行就是计算合并之后的rank值了，而合并之后的rank值应该存在x[]里面，但我们计算的时候又必须用到上一层的rank值，也就是现在x[]
            # 里面放的东西，如果我既要从x[]里面拿，又要向x[]里面放，怎么办？当然是先把x[]的东西放到另外一个数组里面，即y[]。
            x, y = y, x
            p = 1
            x[sa[0]] = 0
            for i in range(1, n):
                if ((y[sa[i - 1]] == y[sa[i]]) and (y[sa[i - 1] + k] == y[sa[i] + k])):
                    x[sa[i]] = p - 1
                else:
                    x[sa[i]] = p
                    p += 1
            # 这里就是用x[]存储计算出的各字符串rank的值了，记得我们前面说过，计算sa[]值的时候如果字符串相同是默认前面的更小的，
            # 但这里计算rank的时候必须将相同的字符串看作有相同的rank，要不然p == n之后就不会再循环啦。
            m = p
            k <<= 1  # k*2

        self.B = ""
        for i in sa:
            self.B += self.T_[(i + self.len) % (self.len + 1)]
        self.S = sa

    def BWT(self):
        self.bz_sort()
        # print(self.B)
        # print(self.S)
        self.O = []  # O(a,i)：在B[0,i]的字符上，a的出现次数。
        a = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '#': 0}
        a[self.B[0]] += 1
        self.O.append(a)
        for i in range(1, len(self.B)):
            a = self.O[-1].copy()
            a[self.B[i]] += 1
            self.O.append(a)
        self.num_a = a['A']
        self.num_c = a['C']
        self.num_g = a['G']
        self.num_t = a['T']

        self.begin_a = 1
        self.begin_c = self.num_a + 1
        self.begin_g = self.num_a + self.num_c + 1
        self.begin_t = self.num_a + self.num_c + self.num_g + 1

    def _select_prefix(self, sub):
        sub = sub[::-1]
        i = len(sub) - 1
        now = sub[i]
        begin = 0
        end = i

        # 最后一个字符
        if now == 'A':
            begin = self.begin_a
            end = self.begin_c - 1
        elif now == 'C':
            begin = self.begin_c
            end = self.begin_g - 1
        elif now == 'G':
            begin = self.begin_g
            end = self.begin_t - 1
        elif now == 'T':
            begin = self.begin_t
            end = self.len

        # 从后往前匹配
        i -= 1
        while i >= 0:
            now = sub[i]
            lbegin = self.B.find(now, begin, end + 1)
            lend = self.B.rfind(now, begin, end + 1)
            if lbegin == -1:
                break

            lnum = self.O[lbegin][now]
            rnum = self.O[lend][now]
            if now == 'A':
                begin = self.begin_a + lnum - 1
                end = self.begin_a + rnum - 1
            elif now == 'C':
                begin = self.begin_c + lnum - 1
                end = self.begin_c + rnum - 1
            elif now == 'G':
                begin = self.begin_g + lnum - 1
                end = self.begin_g + rnum - 1
            elif now == 'T':
                begin = self.begin_t + lnum - 1
                end = self.begin_t + rnum - 1
            i -= 1

        i += 1
        starts = []
        length = len(sub) - i
        for i in self.S[begin:end + 1]:
            starts.append(self.len - i - length)
        return starts, length

    # 得到所有s2与后缀树的部分匹配   s2的index，T的后缀号starts，匹配长度length，
    # [[index, length, [starts]]...]
    def _findCommonStrings(self, s2: str):
        index = 0
        results = []
        while index <= len(s2) - 1:
            starts, length = self._select_prefix(s2[index:])  # 对s2的每一个后缀，#匹配后缀号，匹配长度
            if starts and length > 1:  # 后缀号不为空，且匹配长度超过1，可跳！
                results.append([index, length, starts])
                index += length
            else:
                index += 1  # 否则，相当于按步+1
        return results  # [index, length, [starts]]


    # 从[[index, length, [starts]]...] 的[start]中选出一个后缀号
    def _MultiReg(self, results):
        start = float('Inf')  # 起始正无穷
        end = -float('Inf')  # 结尾负无穷
        length = results[-1][0] + results[-1][1]  # 以最后一项的 index+length_i作为 字符串总长度
        for result in results:  # 找出所有后缀号列表中，最小后缀号做  start起点，最大后缀号+length_i做end终点
            s = min(result[2])  # 最小后缀号
            e = max(result[2]) + result[1]  # 最大后缀号+length_i
            start = s if s < start else start
            end = e if e > end else end
        delete_k = []
        for k, result in enumerate(results):  # 枚举遍历
            if len(result[2]) > 1:  # 如果后缀号列表中后缀号不止一个，选出一个
                rate = result[0] / length  # 每个后缀号：  i = abs(i/(end-start)-index/length)
                temp = np.array(result[2])  # 找出最接近index的？！？？？？？？？？？？？
                temp = (temp / (end - start)) - rate
                temp = list(abs(temp))
                index = temp.index(min(temp))
                result[2] = result[2][index]
            else:  # 如果后缀号列表中后缀号只有一个，不用选，就是这个
                result[2] = result[2][0]
            if result[1] < self.thorehold: delete_k.append(k)
        de_results = [results[i] for i in range(len(results)) if i not in delete_k]
        return de_results


    def _trace_back(self, results: list, p: list) -> list:
        # 依据p，按相同规则， 回溯找出 select_results 对应的索引
        i = p.index(max(p))
        track = [i]
        while i > 0:
            j = i - 1
            if p[i] == results[i][1]: break
            while j >= 0:
                if results[i][2] >= (results[j][2] + results[j][1]) and p[i] == p[j] + results[i][1]:
                    track.append(j)
                    i = j
                    break
                j -= 1
        return track[::-1]

    def _select_CommonStrings(self, s2): #依据动态规划，选出合适的不重叠的同源区段
        results = self._findCommonStrings(s2)  # [[index, length, [starts]]...]
        # return self.choose(results), 0
        select_results = self._MultiReg(results)  # 从[[index, length, [starts]]...] 的[start]中选出一个后缀号
        # return self.choose2(select_results), 0
        m = len(select_results)
        if m == 0: return [], 0
        if m == 1: return select_results, select_results[0][1] / self.len
        p = [i[1] for i in select_results]
        for i in range(1, m):
            for j in range(i):
                if select_results[i][2] >= (select_results[j][2] + select_results[j][1]):
                    p[i] = max(p[i], p[j] + select_results[i][1])
        selected_results = [select_results[i] for i in self._trace_back(select_results, p)]
        return selected_results, max(p) / self.len

    def align(self, s2):
        results, _ = self._select_CommonStrings(s2)  # 无冗余的 [[S_index, length, T_starts]...]
        A, B = [], []  # 比对前分段列表
        align_A, align_B = [], []  # 比对结果，串列表
        s_a, s_b = 0, 0  # 下一段，起始索引
        for r in results:
            A.append(self.T[s_a:r[2]])  # A收集T串 未配对区域，T[s_a:T_starts]
            B.append(s2[s_b:r[0]])  # B收集S串 未配对区域，S[s_b:S_index]
            s_a = r[1] + r[2]  # T串下一段，起始索引
            s_b = r[1] + r[0]  # S串下一段，起始索引
        A.append(self.T[s_a:-1])  # A添加尾余串
        B.append(s2[s_b:])  # B添加尾余串
        for i in range(len(A)):
            _, temp_A, temp_B = PSA_AGP_Kband(A[i], B[i])  # 对A,B 每一段，分别进行PSA_AGP_Kband比对
            align_A.append(temp_A)  # align_A比对结果串列表
            align_B.append(temp_B)  # align_B比对结果串列表
        s_A, s_B = '', ''  # 比对结果串
        for j in range(len(results)):
            s_A += align_A[j]  # s_A收集T串 未配对区域
            s_B += align_B[j]  # s_B收集S串 未配对区域
            s_A += s2[results[j][0]:(results[j][0] + results[j][1])]  # s_A收集T串 配对区域
            s_B += s2[results[j][0]:(results[j][0] + results[j][1])]  # s_B收集S串 配对区域
        s_A += align_A[-1]  # A添加尾余串
        s_B += align_B[-1]  # B添加尾余串
        score = Compute_two(s_A, s_B)  # 求取比对得分
        print("score :", score)
        return score, s_A, s_B  # 返回比对结果，比对得分

def write(label: list, strs:list, filepath:str):
    """
    write sequences to a specific file

    Args:
        strs: sequences
        filepath: file path
    """
    i = 0
    with open(filepath, 'w') as f:
        for s in strs:
            f.write("> " + label[i] + '\n')
            f.write(s)
            f.write("\n")
            i += 1


if __name__ == '__main__':
    label,(A,B) = readfasta('data/SARS2.fasta')[0:2]
    A_trie = PSA_BWT(A) #先对长序列建立bwt
    s2, A_aligned2, B_aligned2 = A_trie.align(B) #再对短序列在bwt中遍历
    write(label,[A_aligned2,B_aligned2],'data/SARS2_align.fasta')


