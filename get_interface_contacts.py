from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from tqdm import tqdm
import time
from Bio.PDB import *
from collections import defaultdict


# 使用示例
pdb_file = "/home/zjlab/score_set/targets/T029.1.pdb"


threshold = 10.0  # 阈值为10埃
output_file = "/home/zjlab/score_set/Target29/contacts.ifc_seq"
f = open(output_file, 'w')

def convert_three_letter_to_one_letter(three_letter_code):
    code_map = {
    'ALA': 'A',
    'ASN': 'N',
    'ASP': 'D',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'LYS': 'K',
    'LEU': 'L',
    'ILE': 'I',
    'PHE': 'F',
    'TRP': 'W',
    'SER': 'S',
    'THR': 'T',
    'MET': 'M',
    'VAL': 'V',
    'CYS': 'C',
    'TYR': 'Y',
    'PRO': 'P',
    'HIS': 'H',
    'ARG': 'R'
}
        
    return code_map.get(three_letter_code, '')

def run(pdb_file, threshold):
    # 创建PDB解析器
    parser = PDBParser()
    
    # 解析PDB文件
    structure = parser.get_structure("protein", pdb_file)

    
    # 获取第一个模型
    model = structure[0]
    
    # 获取蛋白质链
    chain1 = model['D']  # 第一个链
    chain2 = model['F']  # 第二个链
    
    # 用于存储界面残基的列表
    interface_residues1,interface_residues2= [],[]


    
    for residue1 in chain1:
        for residue2 in chain2:
            # 计算残基之间的距离
            distance = residue1['CA'] - residue2['CA']  # 以Cα原子为参考计算距离，可根据需要修改

            # 判断距离是否小于等于阈值
            if distance <= threshold:
                # 如果距离小于等于阈值，则将残基标记为界面残基
                residue1.id = (' ', residue1.id[1], ' ')
                residue2.id = (' ', residue2.id[1], ' ')
                interface_residues1.append(residue1)
                interface_residues2.append(residue2)
                
                first_type=convert_three_letter_to_one_letter(residue1.resname)
                first_id=residue1.id[1]
                

                second_id=convert_three_letter_to_one_letter(residue2.resname)
                second_type=residue2.id[1]
                
                f.write(f"{first_type}{first_id}-{second_type}{second_id}\n")
                
    
   
         
         
    # 用于存储界面残基的列表
    interface_residues = [] 
         
    for residue2 in chain2:
        for residue1 in chain1:
            # 计算残基之间的距离
            distance = residue2['CA'] - residue1['CA']  # 以Cα原子为参考计算距离，可根据需要修改

            # 判断距离是否小于等于阈值
            if distance <= threshold:
                # 如果距离小于等于阈值，则将残基标记为界面残基
                residue2.id = (' ', residue2.id[1], ' ')
                interface_residues.append(residue2) 
                
    # 打印界面残基的信息
    print("Interface residues:")
    for residue in interface_residues:
        print("Chain {}, Residue {}:{}".format('F', convert_three_letter_to_one_letter(residue.resname), residue.id[1]))
    
    print(len(interface_residues))
                
'''
first_id=residue1.get_id()[1]-1
first_type=residue1.get_resname()[:1]

second_id=residue2.get_id()[1]-1
second_type=residue2.get_resname()[:1]
if first_id<=len(chain1) and second_id<=len(chain2):
    f.write(f"{first_type}{first_id}-{second_type}{second_id}\n")
'''
                


print(f'{time.time()}  Start')
# 计算界面接触
run(pdb_file, threshold)


print(f"{time.time()}  calculate_interface_residues Finished!")