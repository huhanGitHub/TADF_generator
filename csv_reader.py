# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: csv_reader.py
# @Author: han
# @Institution: Monash University, Melbourne, Australia
# @E-mail: huhanthu@qq.com, han.hu@monash.edu
# @Site: 
# @Time: 12 18, 2020
# ---
import os

csv_path = 'data\main.csv'
save_dir = 'data\he_smi'

def read_csv(csv_path, save_dir):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    with open(csv_path, 'r', encoding='utf8') as f:
        line_count = 0
        for line in f.readlines():
            line_count += 1
            if line_count > 5:
                file_path = os.path.join(save_dir, str(line.split(',')[0]) + '.smi')
                smi = line.split(',')[1]
                with open(file_path, 'a+', encoding='utf8') as new_f:
                    new_f.write(smi)


if __name__ == '__main__':
    read_csv(csv_path, save_dir)