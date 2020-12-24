# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: excel_data_split.py
# @Author: han
# @Institution: Monash University, Melbourne, Australia
# @E-mail: huhanthu@qq.com, han.hu@monash.edu
# @Site: 
# @Time: 12 23, 2020
# ---

from rdkit.Chem import Draw
import xlrd
import shutil
import os
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit import Chem

smis_dir = 'D:\Machine_Learning\data\smi\smi'
excel_src = 'D:\Machine_Learning\data\data_processed\processed-data.xlsx'
dst = 'D:\Machine_Learning\data\data_processed\EL'

selected_train_data = 'D:\Machine_Learning\data\data_processed\el_selected\\train'
selected_val_data = 'D:\Machine_Learning\data\data_processed\el_selected\\val'
selected_test_data = 'D:\Machine_Learning\data\data_processed\el_selected\\test'

def convert_smis(dir):
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)
        if os.path.isdir(file_path):
            convert_smis(file_path)
        else:
            if file.endswith('smi'):
                with open(file_path) as f:
                    line= f.read()
                    smi=line.split('\t')[0]
                    m = Chem.MolFromSmiles(smi)
                    try:
                        img = Draw.MolToImage(m, options=opts)
                        save_path = file_path + '.png'
                        img.save(save_path)
                    except ValueError:
                        print('Null molecule provided: ' + file)


def excel_data_split(src, dst, smis_dir):
    category_dic = {}
    workbook = xlrd.open_workbook(src)
    sheet1 = workbook.sheet_by_name('Sheet1')
    print(sheet1.name, sheet1.nrows, sheet1.ncols)
    for i in range(1, sheet1.nrows):
        try:
            el = int(sheet1.cell(i, 7).value)
            file = str(sheet1.cell(i, 13).value).strip() + '.smi' + '.png'
            tadf_path = os.path.join(smis_dir, file)

            # dense between 22 and 29
            category = int(el / 10)
            category_dir = os.path.join(dst, str(category-44))
            if not os.path.exists(category_dir):
                os.mkdir(category_dir)
            dst_path = os.path.join(category_dir, file)
            shutil.copy(tadf_path, dst_path)

        except FileNotFoundError as e1:
            print('no file: '+str(e1))
        except ValueError as e2:
            print('value error: ' + str(e2))
        except OSError:
            print('OSError: [Errno 22] Invalid argument: ' + tadf_path)


def test_data_split(dir, ratio, test_dir):
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    test_case_cat = []
    for category in os.listdir(dir):
        cat_dir = os.path.join(dir, category)
        files = []
        for file in os.listdir(cat_dir):
            file_path = os.path.join(cat_dir, file)
            files.append(file_path)
        test_files = files[int(ratio * len(files)):]
        # test_cat_dir = os.path.join(test_dir, category)

        for test_file in test_files:
            file_name = test_file.split('\\')[-1]
            test_file_path = os.path.join(test_dir, file_name)
            test_case = str(file_name.split('.')[0]) + ',' + str(int(category)) + '\n'
            test_case_cat.append(test_case)
            shutil.move(test_file, test_file_path)

    test_case_cat_total = os.path.join(test_dir, 'test_result.txt')
    with open(test_case_cat_total, 'a+', encoding='utf8') as f:
        for case in test_case_cat:
            f.write(case)


def val_data_split(dir, ratio, val_dir):
    if not os.path.exists(val_dir):
        os.mkdir(val_dir)
    for category in os.listdir(dir):
        cat_dir = os.path.join(dir, category)
        files = []
        for file in os.listdir(cat_dir):
            file_path = os.path.join(cat_dir, file)
            files.append(file_path)
        val_files = files[int(ratio * len(files)):]
        val_cat_dir = os.path.join(val_dir, category)
        if not os.path.exists(val_cat_dir):
            os.mkdir(val_cat_dir)
        for val_file in val_files:
            file_name = val_file.split('\\')[-1]
            val_file_path = os.path.join(val_cat_dir, file_name)
            shutil.move(val_file, val_file_path)


if __name__ == '__main__':
    # excel_data_split(excel_src, dst, smis_dir)
    test_data_split(selected_train_data, 0.8, selected_test_data)
    val_data_split(selected_train_data, 0.8, selected_val_data)