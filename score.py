'''
compute the SP score of seqs

author: Juntao Chen
time: 12.27.2021
'''

def spscore(strs:list, d=3, e=1, m=1, mis=-2):
    """
    compute the SP score of seqs

    Args:
        strs: sequences list
           d: open gap penalty
           e: extend gap penalty
           m: match score
         mis: mismatch score
    
    Returns:
        score: sp score
    """
    nums = len(strs)
    score = 0
    for i in range(nums):
        for j in range(i+1, nums):
            score += _spTwo(strs[i], strs[j])
    
    return score


def _spTwo(s1:str, s2:str, d=3, e=1, m=1, mis=-2):
    """
    compute the SP score of two seqs

    Args:
        s1: sequence 1
        s2: sequence 2
        d: open gap penalty
        e: extend gap penalty
        m: match score
        mis: mismatch score
    
    Returns:
        score: spscore
    """
    if len(s1) != len(s2):
        print(s1, s2)
        raise ValueError("seqs not aligned!")
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


