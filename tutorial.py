from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import os

opts = DrawingOptions()
opts.includeAtomNumbers = True

def converter(file_name):
    mols = [ mol for mol in Chem.SDMolSupplier(file_name) ]
    outname = file_name.split(".sdf")[0] + ".smi"
    out_file = open( outname, "w" )
    for mol in mols:
        smi = Chem.MolToSmiles(mol)
        name = mol.GetProp("_Name")
        out_file.write( "{}\t{}\n".format(smi, name ))
    out_file.close()


def read_DA_corpus(a_pool, d_pool):
    a_list = []
    d_list = []
    with open(a_pool) as f:
        for line in f.readlines():
            a_list.append(line.replace('\n', ''))
    with open(d_pool) as f:
        for line in f.readlines():
            d_list.append(line.replace('\n', ''))
    return a_list, d_list


def compose():
    a_pool = 'data\\acceptor_pool.txt'
    d_pool = 'data\donor_pool.txt'
    A1 = Chem.SDMolSupplier('data\D_A\A1.sdf')
    # A2 = Chem.SDMolSupplier('data\D_A\A2.sdf')
    # A3 = Chem.SDMolSupplier('data\D_A\A3.sdf')
    # A4 = Chem.SDMolSupplier('data\D_A\A4.sdf')
    # A5 = Chem.SDMolSupplier('data\D_A\A5.sdf')
    #
    # B1 = Chem.SDMolSupplier('data\D_A\B1.sdf')
    # B2 = Chem.SDMolSupplier('data\D_A\B2.sdf')

    D1 = Chem.SDMolSupplier('data\D_A\D1.sdf')
    # D2 = Chem.SDMolSupplier('data\D_A\D2.sdf')
    # D3 = Chem.SDMolSupplier('data\D_A\D3.sdf')

    a_list, d_list = read_DA_corpus(a_pool, d_pool)

    BRICS.BRICSDecompose('c1ccccc1OCCOC(=O)CC')


    fragms = [Chem.MolFromSmiles(), Chem.MolFromSmiles('c1ccccc1OCCOC(=O)CC')]
    ms = BRICS.BRICSBuild(fragms)
    results = []
    for i in ms:
        results.append(i)
    # prods = [next(ms) for x in range(1)]
    # [prod.UpdatePropertyCache(strict=False) for prod in prods]
        Draw.MolsToGridImage(results, molsPerRow=4, subImgSize=(200, 200))


def test():
    smi = 'C=CC(=O)N1CCC(CC1)C2CCNC3=C(C(=NN23)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=O)N'
    m = Chem.MolFromSmiles(smi)
    frags = (BRICS.BRICSDecompose(m))
    mols = []
    for fsmi in frags:
        mols.append(Chem.MolFromSmiles(fsmi))
    #Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols]).show()

    newms = BRICS.BRICSBuild(mols)
    newms = list(newms)
    print(len(newms))
    new = newms[:5]
    #frags = [Chem.MolFromSmiles(x) for x in newms[:5]]
    Draw.MolsToGridImage(new, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in new]).show()


def test2():
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    #Draw.MolToImage(m, options=opts).show()
    m = Chem.MolFromSmiles('COc1c(Br)cccc1OC')
    print("m1 Smiles:", Chem.MolToSmiles(m))
    Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('OC')
    repsmis = ['F', 'Cl', 'Br', 'O']
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols]).show()


def test3():
    m = Chem.MolFromSmiles('c1(ccc2c(c1)c1c(n2c2ccccc2)ccc(c1)[Ne])[Ne]')
    print("m1 Smiles:", Chem.MolToSmiles(m))
    Draw.MolToImage(m, options=opts).show()
    #m = Chem.AddHs(m)
    acceptor_smi_h = Chem.MolToSmiles(m)
    print("m2 Smiles:", acceptor_smi_h)
    #Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('[Ne]')
    donor = Chem.MolFromSmiles('N(c1ccccc1)(c1ccccc1)[H]')
    #donor = Chem.AddHs(donor)
    Draw.MolToImage(donor, options=opts).show()
    repsmis=[Chem.MolToSmiles(donor)]
    #repsmis = ['c1cc2c(cc1)c1c([nH]2)cccc1', 'c1cc(ccc1)Nc1ccccc1', 'c1cc2c(cc1)Oc1c(N2)cccc1']
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]
    smis = list(set(smis))
    # smis.remove(acceptor_smi_h)
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols]).show()


def test4():
    m = Chem.MolFromSmiles('[H]c1c([H])c(C#N)c(-c2ccc3c(c2)[nH]c2ccccc23)c([H])c1C#N')
    print("m1 Smiles:", Chem.MolToSmiles(m))
    Draw.MolToImage(m, options=opts).show()
    m = Chem.AddHs(m)
    print("m2 Smiles:", Chem.MolToSmiles(m))
    Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('[H]')
    # repsmis = ['F', 'Cl', 'Br', 'O']
    rep = Chem.MolFromSmiles('c1cc2c(cc1)c1c([nH]2)cccc1')
    rep = Chem.AddHs(rep)
    rep = Chem.MolToSmiles(rep)
    repsmis = [rep]
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
    res = AllChem.ReplaceSubstructs(m, patt, rep)
    mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]
    smis = list(set(smis))
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols]).show()


def test5():
    m = Chem.MolFromSmiles('COc1c(Br)cccc1OC')
    Draw.MolToImage(m, options=opts).show()
    patt = Chem.MolFromSmarts('OC')
    repsmis = ['F', 'Cl', 'Br', 'O']
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)
    smis = [Chem.MolToSmiles(mol) for mol in mols]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols]).show()


if __name__ == '__main__':
    test5()