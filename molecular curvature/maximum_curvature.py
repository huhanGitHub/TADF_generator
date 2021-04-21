from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions  # Only needed if modifying defaults
from itertools import permutations
import math

opts = DrawingOptions()
opts.includeAtomNumbers = True
opts.bondLineWidth = 2.8


def maximum_curvature(smi):
    mols = []
    m = Chem.MolFromSmiles(smi)
    # Draw.MolToImage(m, options=opts).save('m3d.png')
    # Draw.MolToImage(m).save('m3d1.png')
    mols.append(m)
    m3d = Chem.AddHs(m)
    AllChem.EmbedMolecule(m3d, randomSeed=1)
    # Draw.MolToImage(m3d, size=(250, 250)).show()
    m3d_without_h = Chem.RemoveHs(m3d)
    mols.append(m3d_without_h)
    bonds = m3d.GetBonds()
    single_bonds_id = []

    # find all single bonds without 'H' as end atom
    for bond in bonds:
        bond_type = bond.GetBondType()
        # print(bond_type)
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # print('begin atom: ' + begin_atom.GetSymbol())
        # print('end atom: ' + end_atom.GetSymbol())
        if str(bond_type) == 'SINGLE':
            if end_atom.GetSymbol() != 'H':
                single_bonds_id.append(bond.GetIdx())

    if len(single_bonds_id) == 0:
        print('no single bond, return')
        return
    frags = Chem.FragmentOnBonds(m, single_bonds_id)
    smis = Chem.MolToSmiles(frags)
    smis = smis.split('.')
    frags_ids = []
    for smi in smis:
        # print(smi)
        frag_ids = []
        mols.append(Chem.MolFromSmiles(smi))
        patt = Chem.MolFromSmarts(smi)
        flag = m.HasSubstructMatch(patt)
        if flag:
            atomids = m.GetSubstructMatches(patt)
            # print("matched atom id:", atomids)
            for atomid in atomids:
                frag_ids.append(atomid)

            if frag_ids not in frags_ids:
                frags_ids.append(frag_ids)
        else:
            print("molecular m do not contain group: ", smi)

    # print(frags_ids) frags_ids: [[(14, 15, 16), (17, 26, 27), (36, 37, 38), (39, 48, 49)], [(4, 5, 6, 28)], [(5, 6,
    # 7, 8, 9, 12, 10, 11), (5, 28, 29, 30, 31, 34, 32, 33)] ...]
    frag_angle = []
    for frag_ids in frags_ids:
        for atomids in frag_ids:
            atom_planes = atom_plane_segmentation(atomids, m3d)
            if atom_planes is None or len(atom_planes) == 0:
                continue
            else:
                plane_angles = []
                for i in range(len(atom_planes)):
                    for j in range(i + 1, len(atom_planes)):
                        plane_angle = angle_of_two_planes(atom_planes[i], atom_planes[j])
                        plane_angles.append(plane_angle)

                max_plane_angle = max(plane_angles)
                frag_angle.append(max_plane_angle)

    print('max angle:', max(frag_angle))
    # img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(400, 400), legends=['' for x in mols], options=opts)
    # img.show()
    return frag_angle


def atom_plane_segmentation(atomids, m3d):
    atom_planes = []
    if len(atomids) <= 3:
        return None
    else:
        for i in range(len(atomids)):
            for j in range(i + 1, len(atomids)):
                for k in range(j + 1, len(atomids)):
                    # determine if plane can generate a valid plane: if three points are all in one line
                    plane = [atomids[i], atomids[j], atomids[k]]
                    plane_coordinate = [list(m3d.GetConformer().GetAtomPosition(ii)) for ii in plane]
                    if plane_coordinate[0][0] == plane_coordinate[1][0] == plane_coordinate[2][0] or \
                            plane_coordinate[0][1] == plane_coordinate[1][1] == plane_coordinate[2][1] or \
                            plane_coordinate[0][2] == plane_coordinate[1][2] == plane_coordinate[2][1]:
                        continue

                    atom_planes.append(plane_coordinate)

    return atom_planes


def angle_of_two_planes(plane1, plane2):
    # obtain the normal vector of each plane
    plane1_normal_vector = get_plane_normal_vector(plane1)
    plane2_normal_vector = get_plane_normal_vector(plane2)

    x1, y1, z1 = plane1_normal_vector
    x2, y2, z2 = plane2_normal_vector

    # Dihedral Angle formula for plane vectors
    cos_theta = ((x1 * x2) + (y1 * y2) + (z1 * z2)) / (
            ((x1 ** 2 + y1 ** 2 + z1 ** 2) ** 0.5) * ((x2 ** 2 + y2 ** 2 + z2 ** 2) ** 0.5))

    # math.degrees: convert radians to angles, math.acos: return the arc cosine radian value of this number
    degree = math.degrees(math.acos(cos_theta))
    # dihedral angle is in [0°,180°], but the angle between two planes is in [0°,90°]
    if degree > 90:
        degree = 180 - degree

    # keep 5 digits
    return str(round(degree, 5))


def get_plane_normal_vector(plane):
    A, B, C = plane
    AB_vector = [B[0] - A[0], B[1] - A[1], B[2] - A[2]]  # vector AB
    AC_vector = [C[0] - A[0], C[1] - A[1], C[2] - A[2]]  # vector AC
    B1, B2, B3 = AB_vector
    C1, C2, C3 = AC_vector
    normal_plane = [B2 * C3 - C2 * B3, B3 * C1 - C3 * B1, B1 * C2 - C1 * B2]
    return normal_plane


if __name__ == '__main__':
    # smi = 'CC(C)OC(=O)C(C)NP(=O)(OCC1C(C(C(O1)N2C=CC(=O)NC2=O)(C)F)O)OC3=CC=CC=C3'
    smi = 'c1cccc(N(c2ccc(cc2)c2nc(C#N)c(nc2c2ccccc2)C#N)c2ccc(cc2)c2nc(C#N)c(nc2c2ccccc2)C#N)c1'
    maximum_curvature(smi)
