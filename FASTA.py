'''
to read/write the sequences data from files

author: Juntao Chen
time:12.27.2021
'''
import os
import re


def writefasta(strs:list, filepath:str):
    """
    write sequences to a specific file

    Args:
        strs: sequences
        filepath: file path
    """
    i = 0
    with open(filepath, 'w') as f:
        for s in strs:
            f.write("> " + str(i) + '\n')
            f.write(s)
            f.write("\n")
            i += 1


def readfasta(filename:str, del_ = False):
    """
    read sequences from a specific file

    Args:
        filename: file path
        del_: True/False, remove gaps in sequence

    Returns:
        labels

        strings
    """
    file_path = os.path.join(os.getcwd(), filename)
    with open(file_path, 'r') as f:
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
        strs.append(temp)
    return labels, strs[1:]
