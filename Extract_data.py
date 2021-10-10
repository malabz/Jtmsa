'''
to extract the dna data from files
author: Juntao Chen
'''
import os
import re
import csv

# to find the gap in string
def find_gap(s):
    loc = 0
    state = 0
    for i in s:
        if i == ' ':
            state = 1
        elif state == 1:
            return loc
        loc += 1


def extract_data(filename='DNA.csv'):
    file_path = os.path.join(os.getcwd(), "data", filename)
    # read the csv data
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        datas = list(reader)
    strs = []
    for data in datas[:10000]:
        loc = find_gap(data[0])
        tmp = data[0][loc:]
        tmp_list = tmp.split(" ")
        tmp = ("").join(tmp_list)
        strs.append(tmp)
    return strs

def ex_nex(filename):
    file_path = os.path.join(os.getcwd(), "data", filename)
    with open(file_path, 'r') as f:
        data = [re.split('\s+', line[:-1]) for line in f]
   
    file_path = os.path.join(os.getcwd(), "data", filename[:-3]+'fasta')
    with open(file_path, 'w') as f:
        for d in data:
            f.write("> " + d[0] + '\n')
            f.write(re.sub('-','',d[1]) + '\n')


def write_fasta(strs, filepath):

    i = 0
    with open(filepath, 'w') as f:
        for s in strs:
            f.write("> " + str(i) + '\n')
            f.write(s)
            f.write("\n")
            i += 1
    return 1


def read_fasta(filename="genome.fasta", del_ = False):
    file_path = os.path.join(os.getcwd(), "data", filename)
    with open(file_path, 'r') as f:
        temp = ""
        strs = []
        for line in f:
            if line.startswith('>'):
                if del_:
                    temp = re.sub('-', '', temp)
                strs.append(temp)
                temp = ""
                continue
            if line.startswith('\n'):
                continue
            temp += line[:-1]
        strs.append(temp)
    return strs[1:]

if __name__ == "__main__":
    names = os.listdir(os.path.join(os.getcwd(), "data"))
    names = [i for i in names if 'nex' in i]
    for i in names: ex_nex(i)
