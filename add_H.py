from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import os

dir = 'data\smi'
sdf_path = 'data\sdf\TADF-A-9.sdf'

def converter(file_name):
    mols = [ mol for mol in Chem.SDMolSupplier(file_name) ]
    outname = file_name.split(".sdf")[0] + ".smi"
    out_file = open( outname, "w" )
    for mol in mols:
        smi = Chem.MolToSmiles(mol)
        name = mol.GetProp("_Name")
        out_file.write( "{}\t{}\n".format(smi, name ))
    out_file.close()

def add_H(file_path):
    print(file_path)
    save_file_path = file_path + '_add_H.smi'
    with open(file_path, 'r', encoding='utf-8') as f, open(save_file_path, 'a+', encoding='utf-8') as save_file:
        lines = f.readlines()
        for line in lines:
            line = line.split('\t')
            smi = line[0]
            name = line[1]
            #print(smi + ' ' + name)
            m = Chem.MolFromSmiles(smi)
            if m == None:
                print('error!!!!!!!!!!!!!! + ' + name)
                continue
            m = Chem.AddHs(m)
            smi_h = Chem.MolToSmiles(m)
            line_add_H = smi_h + '\t' + name
            save_file.write(line_add_H)



if __name__ == '__main__':
    # for root, dirs, files in os.walk(dir):
    #     for file in files:
    #         path = os.path.join(root, file)
    #         add_H(path)

    converter(sdf_path)