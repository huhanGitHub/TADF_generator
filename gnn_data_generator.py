from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

atom_pool = ['C', 'N', 'O', 'S', 'H', 'F', 'Cl', 'Br', 'I', 'Se', 'Te', 'Si', 'P', 'B', 'Sn', 'Ge']
# hybridization S:0, SP:1, SP2:2, SP3:3, SP3D:4, SP3D2:5, UNSPECIFIED:6, OTHER:7
hybridization_pool = ['S', 'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'UNSPECIFIED', 'OTHER']
test = set()


def atom_feature_generator(smi, save_list, tag=0):
    edges = []
    mol = Chem.MolFromSmiles(smi)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    for bond in bonds:
        # 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE, 3.0 for TRIPLE
        bond_type = bond.GetBondTypeAsDouble()
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        atom1 = atom1 if atom1 <= atom2 else atom2
        edge = (atom1, atom2, bond_type)
        edges.append(edge)

    atom_features = []
    for atom in atoms:
        #print(atom.GetIdx())
        symbol = atom.GetSymbol()
        index = 0
        if symbol in atom_pool:
            index = atom_pool.index(symbol)
        # num_degree = atom.GetTotalDegree()
        num_H = atom.GetTotalNumHs()
        num_neighbors = len(atom.GetNeighbors())
        aromaticity = atom.GetIsAromatic()
        aromaticity = 1 if aromaticity is True else 0
        # hybridization S:0, SP:1, SP2:2, SP3:3, SP3D:4, SP3D2:5, UNSPECIFIED:6, OTHER:7
        hybridization = hybridization_pool.index(str(atom.GetHybridization()))
        ring = atom.IsInRing()
        ring = 1 if ring is True else 0
        formal_charge = atom.GetFormalCharge()
        atom_feature = (index, num_H, num_neighbors, aromaticity, hybridization, ring, formal_charge)
        atom_features.append(atom_feature)

    # save features
    features = [tag, atom_features, edges]
    save_list.append(features)


def read_csv(path):
    df = pd.read_csv(path)
    print(df)


def features_generator():
    Chromophores_features_save_path = 'Chromophores_features.csv'
    Solvent_features_save_path = 'Solvent_features.csv'
    Chromophores_save_list = []
    Solvent_save_list = []
    data = pd.read_csv('data/DB for chromophore_Sci_Data_rev02.csv')
    Chromophores = data['Chromophore']
    Solvent = data['Solvent']
    Quantum_yield = data['Quantum yield']

    index = 0
    for chromophore, solvent, qy in zip(Chromophores, Solvent, Quantum_yield):
        # use a!=a to judge nan
        if qy != qy:
            index += 1
            continue
        atom_feature_generator(chromophore, Chromophores_save_list, tag=index)
        atom_feature_generator(solvent, Solvent_save_list, tag=index)
        index += 1

    df1 = pd.DataFrame(Chromophores_save_list, columns=['Tag', 'Atom Features', 'Edges'])
    df1.to_csv(Chromophores_features_save_path, index=False)

    df2 = pd.DataFrame(Solvent_save_list, columns=['Tag', 'Atom Features', 'Edges'])
    df2.to_csv(Solvent_features_save_path, index=False)


if __name__ == '__main__':
    features_generator()
    # read_csv('features.csv')