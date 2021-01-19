# to extract the data from the csv file
# author: Juntao Chen

import os 
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
    # read the filename
    # data_path = os.getcwd() + "/data"
    # # filenames = os.listdir(data_path)
    # file_path = data_path + "/" + filename

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
    
    # write_fasta(strs, os.getcwd() + "\\DNA.fasta")

    return strs

def write_fasta(strs, filepath):

    i = 0
    with open(filepath, 'w') as f:
        for s in strs:
            f.write("> " + str(i) + '\n')
            f.write(s)
            f.write("\n")
            i += 1
    return 1


def read_fasta(filename="genome.fasta"):
    # read the filename
    # data_path = os.getcwd() + "/data"
    # filenames = os.listdir(data_path)
    # file_path = data_path + "/" + filename

    file_path = os.path.join(os.getcwd(), "data", filename)


    with open(file_path, 'r') as f:
        temp = ""
        strs = []
        for line in f:
            if line.startswith('>'):
                strs.append(temp)
                temp = ""
                continue
            if line.startswith('\n'):
                continue
            temp += line[:-1]
        strs.append(temp)
    return strs[1:]

if __name__ == "__main__":
    data = read_fasta()
    print(len(data), len(data[1]), len(data[-1]), len(data[-2]))
    # extract_data()
    # pass