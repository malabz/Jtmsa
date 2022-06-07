'''
centerAlign for DNA sequences, can control the number of sequences in memory

author: Juntao Chen
date: 6.7.2022
'''

import os
import time
import shutil
import argparse

from PSA_STree import PSA_STree

def coast_time(func):
    def fun(*args, **kwargs):
        t = time.perf_counter()
        result = func(*args, **kwargs)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()), end =  "")
        print(f' COST TIME: {time.perf_counter() - t:.3f}s')
        return result
    return fun

# Part 1: 读取文件并将文件分开存储并确定中心序列
@coast_time
def readAndWriteSmallFiles(filePath:str, chunk_size:int, outDir):
    """
    read sequence and seperate the file into files

    Args:
        filePath: file's path
        chunk_size: the number of sequences in buffer list
        outDir: to store the temp file
    Returns:
        centerSeq

        centerSeqIndex: [index of file, index of seq]

        NumOfSequences
    """
    # 放置临时文件的文件夹
    if (not os.path.exists(outDir)):
        os.makedirs(outDir)
    bufferLabels = []
    bufferStrs = []
    # 缓冲区的序列条数
    count = 0
    index = 0
    # 中心序列
    center_Seq = ""
    center_Idx = []
    with open(filePath, 'r') as f_read:
        temp = ""
        for line in f_read:
            # 记录序列标签
            if (line.startswith('>')):
                bufferLabels.append(line[1:-1])
                # 记录序列
                if (temp != ""):
                    bufferStrs.append(temp)
                    # 记录中心序列
                    if (len(temp) > len(center_Seq)):
                        center_Seq = temp
                        center_Idx = [index, count]
                    temp = ""
                    count += 1
                    # 缓冲区满了就写入到磁盘
                    if (count == chunk_size):
                        writeBufferfasta(bufferStrs, bufferLabels, os.path.join(outDir, str(index) + ".fasta"))
                        index += 1
                        # 清空缓冲区
                        count = 0
                        bufferStrs = []
                        # 因为多读取了一个label，所以需要放回
                        bufferLabels = [bufferLabels[count]]

            # 拼接序列
            elif not line.startswith('\n'):
                temp += line[:-1]
        bufferStrs.append(temp)
        # 记录中心序列
        if (len(temp) > len(center_Seq)):
            center_Seq = temp
            center_Idx = [index, count]
        # 最后一组数据落盘
        writeBufferfasta(bufferStrs, bufferLabels,  os.path.join(outDir, str(index) + ".fasta"))
    
    # 返回中心序列和中心序列的idx和序列的总条数
    return center_Seq, center_Idx, index * chunk_size + count


# 将缓冲区数据写入磁盘                     
def writeBufferfasta(strs: list, labels: list, outPath: str, append = False):
    """
    write sequences in buffer area to a specific file

    Args:
        strs: sequences
        labels: labels
        outDir: output dir
        index: index
    """
    i = 0
    writeOrappend = 'a' if append else 'w'
    with open(outPath, writeOrappend) as f:
        for s in strs:
            f.write(">" + labels[i] + '\n')
            for idx in range(0, len(s)//60+1):
                f.write(s[idx*60: (idx+1)*60] + '\n')
            i += 1


# 读取数据
def readfasta(filePath: str):
    """
    read sequences from a specific file

    Args:
        filename: file path

    Returns:
        labels

        strings
    """
    with open(filePath, 'r') as f:
        temp = ""
        strs = []
        labels = []
        for line in f:
            if line.startswith('>'):
                labels.append(line[1:-1])
                strs.append(temp)
                temp = ""
                continue
            if line.startswith('\n'):
                continue
            temp += line[:-1]
        strs.append(temp)
    return labels, strs[1:]


# Part 2: 将小文件一一读取并比对
@coast_time
def AlignToSmallFile(cenSeq: str, numFile: int, outDir: str):
    stree = PSA_STree(cenSeq)
    marks = [0] * (len(cenSeq) + 1)
    for i in range(0, numFile):
        labels, strs = readfasta(os.path.join(outDir, str(i) + ".fasta"))
        strsed_cen = []
        strsed_oth = []
        for string in strs:
            _, _, sA, sB = stree.align(string)
            strsed_cen.append(sA)
            strsed_oth.append(sB)
        # 将比对好的序列写回磁盘
        writeBufferfasta(strsed_oth, labels, os.path.join(outDir, str(i) + ".fasta"))
        writeBufferfasta(strsed_cen, labels, os.path.join(outDir, str(i) + ".cen.fasta"))
        # 得到gap矩阵
        marks = getMarks(strsed_cen, marks)
    return marks;


def getMarks(strsed: list, marks: list):
    for string in strsed:
        i = 0
        counter = 0
        for c in string:
            if c == '-':
                counter += 1
            else:
                marks[i] = max(marks[i], counter)
                counter = 0
                i += 1
        marks[i] = max(marks[i], counter)
    return marks


# Part 3: 插入gap，完成比对
@coast_time
def insertGapsToAll(marks: list, numFile: int, outDir: str, outPath: str):
    for i in range(0, numFile):
        labels, strs = readfasta(os.path.join(outDir, str(i) + ".fasta"))
        _, strs_cen = readfasta(os.path.join(outDir, str(i) + ".cen.fasta"))
        newStrs = []
        for idx in range(len(strs)):
            newmark = [0] * (len(strs_cen[idx]) + 1)
            pi = 0
            pj = 0
            total = 0
            for c in strs_cen[idx]:
                if c == '-':
                    total += 1
                else:
                    newmark[pi] = marks[pj] - total;
                    pi += 1
                    pj += 1
                    while total != 0:
                        pi += 1
                        total -= 1
            newmark[pi] = marks[pj] - total;
            newStrs.append(insertGap(newmark, strs[idx]))
            # print(len(newStrs[len(newStrs) - 1]), end = ", ")
        writeBufferfasta(newStrs, labels, outPath, True)


def insertGap(mark: list, str: str):
    res = ""
    length = len(mark)
    for i in range(length):
        res += ("-" * mark[i])
        if i < length - 1:
            res += str[i]
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, default='', help="The input Fasta file's path.")
    parser.add_argument('--output', type=str, required=True, default='', help="The output Fasta file's path.")
    parser.add_argument('--nums', type=int, default=100, help='Number of sequences resident in memory, default=100.')
    args = parser.parse_args()

    filePath = args.input
    outPath = args.output
    chunk_size = args.nums
    outDir = os.path.join(os.path.dirname(filePath), ".tmp")

    print("-----------Args List-----------")
    print("  Input Path: " + filePath)
    print(" Output Path: " + outPath)
    print("    tmp Path: " + outDir)
    print("  Chunk Size: " + str(chunk_size))
    print("-----------Args List-----------")
    
    try:
        # 删除可能存在的tmp文件
        if os.path.exists(outDir):
            shutil.rmtree(outDir)
        # 删除已经存在输出文件
        if os.path.exists(outPath):
            os.remove(outPath)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()), end =  "")
        print(" STEP 1: SPLIT the Sequences")
        seq, idx, num = readAndWriteSmallFiles(filePath, chunk_size, outDir)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()), end =  "")
        print(" STEP 2: ALIGN the Sequences")
        marks = AlignToSmallFile(seq, num // chunk_size + 1, outDir)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()), end =  "")
        print(" STEP 3: COLLECT the Sequences")
        insertGapsToAll(marks, num // chunk_size + 1, outDir, outPath)
        print(" DONE!")
    finally:
        if os.path.exists(outDir):
            shutil.rmtree(outDir)