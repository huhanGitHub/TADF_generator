# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: run_obabel.py.py
# @Author: han
# @Institution: Monash University, Melbourne, Australia
# @E-mail: huhanthu@qq.com, han.hu@monash.edu
# @Site: 
# @Time: 12 18, 2020
# ---

import os
src = 'data\obabel'
dst = 'pic/'

def run_obabel(src, dst):
    if os.path.isfile(src):
        with open(src, 'r', encoding='utf8') as f:
            smis = f.readlines()
            for i in range(len(smis)):
                smis[i] = smis[i].replace('\n', '')
                dir = str(i % 5000)
                output_dir = os.path.join(dst, dir)
                if not os.path.exists(output_dir):
                    os.mkdir(output_dir)
                output_path = os.path.join(output_dir, str(i)+'.png')
                cmd = 'obabel ' + str(smis[i]) + ' -O ' + output_path
                os.system(cmd)
    else:
        if os.path.isdir(src):
            for file in os.listdir(src):
                run_obabel(file)

if __name__=='__main__':
    run_obabel(src, dst)