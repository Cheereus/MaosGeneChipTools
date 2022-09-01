from tqdm import trange
import numpy as np
import time
import sys


def read_file_nd(filePath, head=0):
    print('Reading...', filePath)
    time.sleep(0.5)

    f = open(filePath, encoding='utf-8')
    lines = f.readlines()
    f.close()
    data = []
    for i in trange(head, len(lines)):
        line = lines[i]
        if '#' in line:
            continue
        if '\n' in line:
            line = line.replace('\n', '')
            # print('回车')
        line = list(map(str, line.split()))
        data.append(line)
    data = np.array(data)
    time.sleep(0.5)
    print(filePath, data.shape)
    return data


def database_initialize(data_list):
    data_dict = {}
    data_duplicate = []
    data_duplicate_dict = {}
    duplicate_list = []
    duplicate_ID_dict = {}
    for data_name in data_list:
        data = read_file_nd(data_name + '.map')
        for line in data:
            line = line.tolist()
            CHR, POS, ID = line[0], line[3], line[1]
            if CHR == '0' or POS == '0':
                continue
            key_name = ':'.join([CHR, POS])
            if key_name in data_dict:
                if ID not in duplicate_ID_dict:
                    duplicate_list.append([ID, data_dict[key_name][1]])
                duplicate_ID_dict[ID] = True

                if data_dict[key_name][3] not in data_duplicate_dict:
                    line2 = data_dict[key_name]
                    if data_name == line2[-1]:
                        data_duplicate.append(line + [data_name])
                        data_duplicate.append(line2)
                    data_duplicate_dict[ID] = True
            else:
                data_dict[key_name] = line + [data_name]

    keys = list(data_dict.keys())
    output1 = open('database.txt', 'w', encoding='utf-8')
    for i in trange(len(keys)):
        line = data_dict[keys[i]]
        output1.writelines('\t'.join(line) + '\n')
    output1.close()

    output2 = open('duplicate.txt', 'w', encoding='utf-8')
    for i in trange(len(data_duplicate)):
        line = data_duplicate[i]
        output2.writelines('\t'.join(line) + '\n')
    output2.close()

    output3 = open('duplicate_2.txt', 'w', encoding='utf-8')
    for i in trange(len(duplicate_list)):
        line = duplicate_list[i]
        output3.writelines('\t'.join(line) + '\n')
    output3.close()


if __name__ == '__main__':
    cp_l = sys.argv[1:]
    database_initialize(cp_l)









