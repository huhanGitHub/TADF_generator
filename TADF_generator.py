from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import os

def read_donor_acceptor_corpus(a_pool, d_pool):
    a_list = []
    d_list = []
    with open(a_pool) as f:
        for line in f.readlines():
            a_list.append(line.replace('\n', ''))
    with open(d_pool) as f:
        for line in f.readlines():
            d_list.append(line.replace('\n', ''))
    return a_list, d_list


def concat_donor_acceptor(acceptor_smi, donors_smi, save_dir, count):
    m = Chem.MolFromSmiles(acceptor_smi)
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    #print("m1 Smiles:", Chem.MolToSmiles(m))
    #Draw.MolToImage(m, options=opts).show()
    m = Chem.AddHs(m)
    acceptor_smi_h = Chem.MolToSmiles(m)
    #print("m2 Smiles:", acceptor_smi_h)
    #Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('[H]')
    repsmis = donors_smi
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]

    smis = list(set(smis))
    smis.remove(acceptor_smi_h)

    for smi in smis:
        count = concat_donor_donor_acceptor(smi, donors_smi, save_dir, count)

    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    count = count + len(mols)
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols])
    #img.show()
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    save_path = os.path.join(save_dir, acceptor_smi + '.jpg')
    img.save(save_path)
    return count

def concat_donor_donor_acceptor(acceptor_donor_smi, donors_smi, save_dir, count):
    m = Chem.MolFromSmiles(acceptor_donor_smi)
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    #print("m1 Smiles:", Chem.MolToSmiles(m))
    #Draw.MolToImage(m, options=opts).show()
    m = Chem.AddHs(m)
    acceptor_smi_h = Chem.MolToSmiles(m)
    #print("m2 Smiles:", acceptor_smi_h)
    #Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('[H]')
    repsmis = donors_smi
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]

    smis = list(set(smis))
    smis.remove(acceptor_smi_h)

    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    count = count + len(mols)
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols])
    #img.show()
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    save_path = os.path.join(save_dir, acceptor_donor_smi + '.jpg')
    img.save(save_path)
    return count

if __name__ == '__main__':
    save_dir = 'TADF'
    a_pool = 'data\\acceptor_pool.txt'
    d_pool = 'data\donor_pool.txt'
    count = 0
    a_list, d_list = read_donor_acceptor_corpus(a_pool, d_pool)
    for acceptor in a_list:
        count = concat_donor_acceptor(acceptor, d_list, save_dir, count)
    print(count)


