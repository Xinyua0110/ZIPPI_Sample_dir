from Bio import AlignIO
import os
import csv
from collections import Counter
import math
import numpy as np




def caculated_DI(msa_file1,msa_file2,ifc_seq_file):
    # 读取第一个.msa文件

    alignment1 = AlignIO.read(msa_file1, "fasta")

    # 读取第二个.msa文件
    alignment2 = AlignIO.read(msa_file2, "fasta")

    
    # 读取contacts位点信息
    with open(ifc_seq_file, "r") as f:
        lines = f.readlines()
    contacts = []
    for line in lines:
        aa1, aa2 = line.strip().split("-")
        pos1 = int(aa1[1:])
        pos2 = int(aa2[1:])
        contacts.append((pos1, pos2))

    A_list,B_list=[],[]
    a_list=[]
    b_list=[]
    for i in range(len(contacts)):
        aa_list = []
        for record in alignment1:
            index=contacts[i][0]
            aa = record.seq[index-1]  # 第5个位点的索引是4，因为序列索引从0开始
            aa_list.append(aa)
        A_list.append(aa_list)
    for j in range(len(contacts)):       
        bb_list = []
        for record in alignment2:
            index=contacts[j][1]
            bb = record.seq[index-1]  # 第5个位点的索引是4，因为序列索引从0开始
            bb_list.append(bb)
        B_list.append(bb_list)
        
    # 遍历a_list和b_list中的所有列表
    DI_list=[]
    for i in range(len(A_list)):
        # 从a_list和b_list中取出第i个list
        list_a = A_list[i]
        list_b = B_list[i]
        
        # 将list_a和list_b合并成一个列表
        combined_list = list_a + list_b

        # 计算氨基酸频率
        aa_count = Counter(combined_list)

        
        # 取两个列表中相同位置上的元素，得到氨基酸对列表
        aa_pairs = [(aa_a, aa_b) for aa_a, aa_b in zip(list_a, list_b) if aa_a == aa_b]

        # 计算氨基酸对的频率
        aa_pair_count = Counter(aa_pairs)
        
        # 计算 DI(a,b)
        DI = 0.0
        for aa_a in set(list_a):
            for aa_b in set(list_b):
                p_aa = (aa_count[aa_a] ) / float(len(combined_list))
                p_bb = (aa_count[aa_b] ) / float(len(combined_list))
            if (aa_a, aa_b) in aa_pair_count:
                p_ab = aa_pair_count[(aa_a, aa_b)] / float(len(aa_pairs))
            else:
                p_ab = 1e-9
            if p_aa <= 0 or p_bb <= 0 or p_ab <= 0:
                DI += 0.0
            else:
                DI += p_ab * np.log(p_ab / (p_aa * p_bb))
        DI_list.append(DI)
        
        
        
            
    return DI_list



def caculated_DI_APC(msa_file1, msa_file2, ifc_seq_file):
    # 计算DI列表
    DI_list = caculated_DI(msa_file1, msa_file2, ifc_seq_file)

    # 读取contacts位点信息
    with open(ifc_seq_file, "r") as f:
        lines = f.readlines()
    contacts = []
    for line in lines:
        aa1, aa2 = line.strip().split("-")
        pos1 = int(aa1[1:])
        pos2 = int(aa2[1:])
        contacts.append((pos1, pos2))

    # 计算DI_avg
    DI_avg = np.mean(DI_list)

    # 计算APC修正值
    APC_list = []
    for i in range(len(contacts)):
        a, b = contacts[i]
        connected_pairs_a = [pair for pair in contacts if a in pair]
        connected_pairs_b = [pair for pair in contacts if b in pair]
        DI_values_a = [DI_list[contacts.index(pair)] for pair in connected_pairs_a]
        DI_values_b = [DI_list[contacts.index(pair)] for pair in connected_pairs_b]
        DI_a_avg = np.mean(DI_values_a)
        DI_b_avg = np.mean(DI_values_b)
        APC = DI_a_avg * DI_b_avg / DI_avg
        APC_list.append(APC)

    

    return APC_list


def read_files_from_csv(csv_file):
    big_DI_APC=[]
    fold_path1="/home/zjlab/Sample_dir/IFC/"
    fold_path2="/home/zjlab/Sample_dir/MSA/"
    # 读取csv文件
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 构造文件路径
            pdb_id = row['PDBID']
            uni1 = row['Uni1']
            uni2 = row['Uni2']
            ifc_seq_file = os.path.join(fold_path1, pdb_id + '_AB.ifc_seq')
            msa_file1 = os.path.join(fold_path2, uni1+'.msa')
            msa_file2 = os.path.join(fold_path2, uni2+'.msa')
            DI=caculated_DI(msa_file1,msa_file2,ifc_seq_file)
            big_DI_APC.append(DI)
            
            

            

            # 返回结果
    return big_DI_APC




hete_DI_APC=read_files_from_csv("/home/zjlab/Sample_dir/bacteria_heteroPDBs_50sample.csv")
homo_DI_APC=read_files_from_csv("/home/zjlab/Sample_dir/bacteria_homoPDBs_50sample.csv")


# 将 heterodimer DI 值写入文件
with open("hete_DI_APC.txt", "w") as f:
    for mi_list in hete_DI_APC:
        for mi in mi_list:
            f.write(str(mi) + " ")
        f.write("\n")
        
# 将 homo DI 值写入文件
with open("homo_DI_APC.txt", "w") as f:
    for mi_list in homo_DI_APC:
        for mi in mi_list:
            f.write(str(mi) + " ")
        f.write("\n")