from bs4 import BeautifulSoup
import lxml
import requests
import os
import rarfile
import re
import eventlet  # 导入eventlet模块
eventlet.monkey_patch()  # 添加猴子补丁


def info_extract(dir, save_dir):
    ENERGY_OSCILLATOR_STRENGTH = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if not str(file).endswith('.rar'):
                file_path = os.path.join(root, file)
                with open(file_path) as f:
                    lines = f.readlines()
                    atoms = []
                    for i in range(len(lines)):
                        if lines[i].startswith(' INPUT CARD>C1  '):
                            ii = i
                            while not lines[ii+1].startswith(' INPUT CARD> $END  '):
                                line = lines[ii+1].replace(' INPUT CARD>', '').strip()
                                atoms.append(line)
                                ii = ii + 1
                                if ii + 1 >= len(lines):
                                    break
                        if i + 1 >= len(lines):
                            break
                        if lines[i].startswith(' STATE #   1  ENERGY'):
                            ENERGY = lines[i].replace(' STATE #   1  ENERGY =    ', '').replace('EV', '').strip()
                            i = i + 1
                            OSCILLATOR_STRENGTH = lines[i].replace(' OSCILLATOR STRENGTH =    ', '').strip()
                            i = i + 1
                            LAMBDA_DIAGNOSTIC = lines[i].replace(' LAMBDA DIAGNOSTIC   =    ', '').replace(' (RYDBERG/CHARGE TRANSFER CHARACTER)', '').strip()
                            i = i + 1
                            ENERGY_OSCILLATOR_STRENGTH.append([file, ENERGY, OSCILLATOR_STRENGTH, LAMBDA_DIAGNOSTIC])

                    atoms_path = os.path.join(save_dir, file)
                    with open(atoms_path, 'a+', encoding='utf8') as f:
                        for atom in atoms:
                            atom = re.sub(' +', ' ', atom)
                            atom = atom.split(' ')
                            try:
                                atom.pop(1)
                            except IndexError:
                                print('IndexError: ' + str(atom))
                                continue
                            atom = ' '.join(atom)
                            f.write(atom + '\n')

    save_csv_path = os.path.join(save_dir, 'all.csv')
    with open(save_csv_path, 'a+', encoding='utf8') as f:
        line = ['atom_id', 'ENERGY (EV)', 'OSCILLATOR_STRENGTH', 'LAMBDA_DIAGNOSTIC']
        f.write(','.join(line) + '\n')

        for i in ENERGY_OSCILLATOR_STRENGTH:
            f.write(','.join(i) + '\n')


def rar_decompress(dir):
    for file in os.listdir(dir):
        if file.endswith('.rar'):
            file_path = os.path.join(dir, file)
            rf = rarfile.RarFile(file_path)
            rf.extractall(file)


def urldownload(url, filename=None):
    """
    :param url
    :param filename: ./test.xls
    :return:
    """
    if os.path.exists(filename) and os.path.getsize(filename) != 0:
        print('file has downloaded: ' + filename)
        return

    with eventlet.Timeout(20, False):
        down_res = requests.get(url)
        with open(filename, 'wb') as file:
            file.write(down_res.content)


def url_crawler(url, save_path):
    response = requests.get(url)
    html = response.text
    soup = BeautifulSoup(html, 'lxml')
    pattern = 'Compound_'
    urls = soup.find_all('a')
    with open(save_path, 'a+', encoding='utf8') as f:
        for l in urls:
            href = l.attrs['href']
            if pattern in href:
                full_href = 'http://pubchemqc.riken.jp/' + href
                f.write(full_href + '\n')


def crawler(url, save_dir):
    response = requests.get(url)
    html = response.text
    soup = BeautifulSoup(html, 'lxml')
    pattern = '.td-b3lyp_6-31g+(d).log.xz'
    urls = soup.find_all('a')
    for l in urls:
        href = l.attrs['href']
        #print(type(href))
        # print(href)
        #print(('.b3lyp_6-31g(d).log.xz' in str(href)))
        if pattern in href:
            # urls_list.append(href)
            text = l.get_text()
            file_name = text.split('.')[0]
            save_path = os.path.join(save_dir, file_name + '.zip')
            full_href = 'http://pubchemqc.riken.jp/' + href
            urldownload(full_href, save_path)
    #print(urls_list)


def main_crawler_urls(urls_path, save_dir):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    with open(urls_path, encoding='utf8') as f:
        lines = f.readlines()
        index = 0
        for url in lines:
            index += 1
            print('index: ' + str(index))
            if index > 500:
                print('500, stop')
                return
            url = url.replace('\n', '')
            name = url.replace('http://pubchemqc.riken.jp/', '').replace('.html', '')
            sub_save_dir = os.path.join(save_dir, name)
            if not os.path.exists(sub_save_dir) and not os.path.isdir(sub_save_dir):
                os.mkdir(sub_save_dir)
            crawler(url, sub_save_dir)


if __name__ == '__main__':
    url = 'http://pubchemqc.riken.jp/Compound_000000001_000025000.html'
    save_dir = r'D:\projects\TADF_generator\data\crawler_data'
    info_dir = r'D:\projects\TADF_generator\data\info_dir'
    #crawler(url, save_dir)
    #rar_decompress(save_dir)
    #info_extract(save_dir, info_dir)

    # url_crawler('http://pubchemqc.riken.jp/', r'D:\projects\TADF_generator\data\jp_urls.txt')

    main_crawler_urls(r'D:\projects\TADF_generator\data\jp_urls.txt', r'D:\projects\TADF_generator\data\jp_data')