# encoding: utf8
import os


def read_csv(path, tag, save_dir):
    smis = []
    results = []
    with open(path, encoding='utf8') as f:
        lines = f.readlines()
        for line in lines:
            index, smile = line.split(',')
            smile = str(smile).replace('\n', '')
            if smile not in smis:
                smis.append(smile)
                results.append([index, smile])

    for i in results:
        file_name = tag + '_' + str(i[0]) + '.smi'
        file_path = os.path.join(save_dir, file_name)
        with open(file_path, 'a+', encoding='utf8') as f:
            f.write(i[1] + '\t' + file_name)


if __name__ == '__main__':
    path = r'D:\projects\TADF_generator\data\Alldata_SMILES_v0.1-1.csv'
    tag = 'A'
    save_dir = r'D:\projects\TADF_generator\data\A'
    read_csv(path, tag, save_dir)