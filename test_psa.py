import os
import datetime
from Extract_data import read_fasta
from PSA_Kband import PSA_AGP_Kband
from PSA_suffixtrie import Suffix_Trie
from PSA_LCS import Nakatsu

if __name__ == "__main__":

    name = ['SARS','dog_eye','Homo_sapiens','SCML4','genome']
    for i in range(len(name)):
        # i = 1
        strs = read_fasta(filename=name[i]+'.fasta')[0:2]
        A = strs[0]
        B = strs[1]
        print(name[i],":",len(A),len(B))

        # 1
        begin_time1 = datetime.datetime.now()
        s1, A_aligned1, B_aligned1 = PSA_AGP_Kband(A, B, d = 3, e = 1)
        end_time1 = datetime.datetime.now()
        print('time_kband:', (end_time1-begin_time1).seconds)

        # 2
        begin_time2 = datetime.datetime.now()
        if len(A) >= len(B):
            A_trie = Suffix_Trie(A)
            s2, A_aligned2, B_aligned2 = A_trie.align(B)
        else:
            B_trie = Suffix_Trie(B)
            s2, A_aligned2, B_aligned2 = B_trie.align(A)
        end_time2 = datetime.datetime.now()
        print('time_suffix:', (end_time2-begin_time2).seconds)

        # 3
        begin_time3 = datetime.datetime.now()
        Na_lcs = Nakatsu(A, B)
        s3, A_aligned3, B_aligned3 = Na_lcs.align()
        end_time3 = datetime.datetime.now()
        print('time_LCS:', (end_time3-begin_time3).seconds)
        print('score:', s1, s2, s3)

        file_path = os.path.join(os.getcwd(), "data", name[i] + "_align.fasta")
        with open(file_path, 'w') as f:
            f.write('\n' + '> 0 ' + 'Kband' + ' length: '+ str(len(A)) + '\n')
            f.write(A_aligned1)
            f.write('\n' + '> ' + str(i) + ' length: '+ str(len(B)) + '\n')
            f.write(B_aligned1)

            f.write('\n' + '> 0 ' + 'Suffix' + ' length: '+ str(len(A)) + '\n')
            f.write(A_aligned2)
            f.write('\n' + '> ' + str(i) + ' length: '+ str(len(B)) + '\n')
            f.write(B_aligned2)
            f.write('\n')

            f.write('\n' + '> 0 ' + 'LCS' + ' length: '+ str(len(A)) + '\n')
            f.write(A_aligned3)
            f.write('\n' + '> ' + str(i) + ' length: '+ str(len(B)) + '\n')
            f.write(B_aligned3)
            f.write('\n')

            f.write('\n' + '> ' + 'SP score of Kband:' + str(s1))
            f.write('\n' + '> ' + 'Time of Kband:' + str((end_time1-begin_time1).seconds))

            f.write('\n' + '> ' + 'SP score of Suffix:' + str(s2))
            f.write('\n' + '> ' + 'Time of Suffix:' + str((end_time2-begin_time2).seconds))

            f.write('\n' + '> ' + 'SP score of LCS:' + str(s3))
            f.write('\n' + '> ' + 'Time of LCS:' + str((end_time3-begin_time3).seconds))