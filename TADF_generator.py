from rdkit import Chem
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

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    # for i in range(len(mols)):
    #     img = Draw.MolToImage(mols[i], options=opts)
    #     save_path = os.path.join(save_dir, acceptor_donor_smi + '_' + str(i) + '.jpg')
    #
    #     # save combined TADF
    #     img.save(save_path)
    return count


def concat_donor_acceptor(acceptor_smi, donors_smi, save_dir, count):
    m = Chem.MolFromSmiles(acceptor_smi)

    # add 'H' to SMILE
    m = Chem.AddHs(m)
    acceptor_smi_h = Chem.MolToSmiles(m)

    # set '[H]' as the connector of donors and acceptors
    patt = Chem.MolFromSmarts('[H]')

    # set 'donors_smi' as the candidate donors to replace the connector
    repsmis = donors_smi
    mols = []
    mols.append(m)
    for r in repsmis:
        rep = Chem.MolFromSmarts(r)

        # replace the connector '[H]' with the candidate donors 'donors_smi'
        res = AllChem.ReplaceSubstructs(m, patt, rep)
        mols.extend(res)

    smis = [Chem.MolToSmiles(mol) for mol in mols]

    # remove repeated smiles and acceptor
    smis = list(set(smis))
    smis.remove(acceptor_smi_h)
    count = count + len(smis)


    # set combined acceptors and donors as the new acceptor to iteratively generate TADF once
    for smi in smis:
        count = concat_donor_donor_acceptor(smi, donors_smi, save_dir, count)

    # mols = [Chem.MolFromSmiles(smi) for smi in smis]
    # opts = DrawingOptions()
    # opts.includeAtomNumbers = True
    # #img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=['' for x in mols])
    # if not os.path.exists(save_dir):
    #     os.mkdir(save_dir)
    # for i in range(len(mols)):
    #     img = Draw.MolToImage(mols[i], options=opts)
    #     save_path = os.path.join(save_dir, acceptor_smi + '_' + str(i) + '.jpg')
    #
    #     # save combined TADF
    #     img.save(save_path)
    return count


def show(smi, save_dir):
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    acceptor = Chem.MolFromSmiles(smi)
    img = Draw.MolToImage(acceptor, options=opts)
    save_path = os.path.join(save_dir, smi + '.jpg')
    img.save(save_path)


def tadf_generator(acceptor_smis, donors_smis, save_dir):
    all_smis = []
    for acceptor_smi in acceptor_smis:
        m = Chem.MolFromSmiles(acceptor_smi)

        # add 'H' to SMILE
        m = Chem.AddHs(m)
        acceptor_smi_h = Chem.MolToSmiles(m)

        # set '[H]' as the connector of donors and acceptors
        patt = Chem.MolFromSmarts('[H]')

        # set 'donors_smi' as the candidate donors to replace the connector
        repsmis = donors_smis
        mols = []
        mols.append(m)
        for r in repsmis:
            rep = Chem.MolFromSmarts(r)

            # replace the connector '[H]' with the candidate donors 'donors_smi'
            res = AllChem.ReplaceSubstructs(m, patt, rep)
            mols.extend(res)

        smis = [Chem.MolToSmiles(mol) for mol in mols]

        # remove repeated smiles and acceptor
        smis = list(set(smis))
        smis.remove(acceptor_smi_h)
        all_smis.extend(smis)
    return all_smis



if __name__ == '__main__':
    save_dir = r'D:\projects\TADF_generator\smis'
    a_pool = r'data\acceptor_pool.txt'
    d_pool = r'data\donor_pool.txt'
    count = 0
    a_list, d_list = read_donor_acceptor_corpus(a_pool, d_pool)

    index = 0
    while count < 5000000:
        smis = tadf_generator(a_list, d_list, save_dir)
        count += len(smis)
        print(count)
        save_path = os.path.join(save_dir, str(index) + '.txt')
        with open(save_path, 'a+', encoding='utf8') as f:
            for smi in smis:
                f.write(smi + '\n')
            index += 1
        if count > 1000000:
            smis = smis[:120000]
        a_list = smis


    # for acceptor in a_list:
    #     show(acceptor, 'data\\patent')
    # for donor in d_list:
    #     show(donor, 'data\\patent')
