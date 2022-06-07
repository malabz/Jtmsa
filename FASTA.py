'''
to read/write the sequences data from files

author: Juntao Chen
time:12.27.2021
'''
import os
import re


def writefasta(strs: list, labels, filepath: str):
    """
    write sequences to a specific file

    Args:
        strs: sequences
        filepath: file path
    """
    i = 0
    with open(filepath, 'w') as f:
        for s in strs:
            f.write(">" + labels[i] + '\n')
            for idx in range(0, len(s)//60+1):
                f.write(s[idx*60: (idx+1)*60] + '\n')
            i += 1


def readfasta(filename: str, del_=False):
    """
    read sequences from a specific file

    Args:
        filename: file path
        del_: True/False, remove gaps in sequence

    Returns:
        labels

        strings
    """
    with open(filename, 'r') as f:
        temp = ""
        strs = []
        labels = []
        for line in f:
            if line.startswith('>'):
                labels.append(line[1:-1])
                if del_:
                    temp = re.sub('-', '', temp)
                strs.append(temp)
                temp = ""
                continue
            if line.startswith('\n'):
                continue
            temp += line[:-1]
        if del_:
            temp = re.sub('-', '', temp)
        strs.append(temp)
    return labels, strs[1:]
