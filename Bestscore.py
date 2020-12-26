# Computing sim(s, t)
# Algorithm BestScore
# input: sequence s and t
# output: vector a
# author: Juntao Chen


# gap
g = -1
# input sequences
t = 'TTGCCATT'
s = 'CCAATTTT'

# aligned sequence
Align_s = ['']*len(s)

# mathch 1 mismatch -1
def p(i,j):
    return int(s[i]==t[j]) - int(s[i]!=t[j])

def BestScore(s, t):
    m = len(s)
    n = len(t)
    a = [0] * (n+1)
    
    for j in range(0, n+1):
        a[j] = j * g

    for i in range(1, m+1):
        old = a[0]
        a[0] = i*g
        for j in range(1, n+1):
            temp = a[j]
            a[j] = max(a[j]+g, old+int(s[i-1]==t[j-1]) - int(s[i-1]!=t[j-1]), a[j-1]+g)
            old = temp
    return a

# Divide and conquer strategy
def Align(s, t, a, b, c, d):

    s_align = list(s)
    t_align = list(t)

    if s[a:b] == '' or t[c:d] == '':
        if s[a:b]:
            Align_s[a:b] = ['_']*len(s[a:b])
    else:
        i = (a+b)//2
        # dimension: d-c
        pref_sim = BestScore(s[a:i], t[c:d])
        
        t_rev = t[c:d][::-1]
        s_rev = s[i+1:b][::-1]
        suff_sim = BestScore(s_rev, t_rev)
        suff_sim = suff_sim[::-1]

        posmax = c-1
        typemax = '_'
        Vmax = pref_sim[0] + g + suff_sim[0]

        # 遍历d-c遍找到与i最匹配的j
        for j in range(c, d):
            if pref_sim[j-c] + p(i,j) + suff_sim[j+1-c] > Vmax:
                posmax = j
                typemax = 'N'
                Vmax = pref_sim[j-c] + p(i,j) + suff_sim[j-c+1]

            if pref_sim[j-c+1] + g + suff_sim[j-c+1] > Vmax:
                posmax = j
                typemax = '_'
                Vmax = pref_sim[j-c] + g + suff_sim[j-c]

        if typemax == '_':
            Align_s[i] = '_' 
            Align(s, t, a, i, c, posmax+1)
            Align(s, t, i+1, b, posmax+1, d)

        else:
            Align_s[i] = posmax
            Align(s, t, a, i, c, posmax)
            Align(s, t, i+1, b, posmax+1, d)

        if '' not in Align_s:
            j = 0
            k_t = 0
            k_s = 0
            for i in Align_s:
                if i == '_':
                    t_align.insert(j+1+k_t, '_')
                    k_t += 1
                elif i-j > 1:
                    if j not in Align_s:
                        j -= 1
                    for ii in range(i-j-1):
                        s_align.insert(ii+Align_s.index(i)+k_s,'_')
                    k_s += i-j-1
                    j = i
                else:
                    j = i

    return s_align, t_align

if __name__ == "__main__":

    Score = BestScore(s, t)
    s_a, t_a = Align(s,t,0,len(s),0,len(t))
    print('Best_score : ', Score[-1])
    print(' '.join(s_a))
    print(' '.join(t_a))