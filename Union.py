# coding=utf-8
# samtools faidx reference/pig.fa chr1:1111-1111
import subprocess
from tqdm import trange
import time
import datetime
import numpy as np
import sys
from functools import wraps
import os


class Union:

    # chips: 芯片名称列表
    def __init__(self, chips):

        # 在此配置各软件包的调用路径
        self.PLINK_PATH = 'plink'
        self.SAMTOOLS_PATH = 'samtools'
        self.BEAGLE_PATH = '/home/liujf/WORKSPACE/maorh/software/beagle.18May20.d20.jar'
        self.PIG_FA_PATH = '/home/liujf/WORKSPACE/maorh/software/pig.fa'

        # 获取任务开始的时间，并以此时间创建日志文件
        
        self.TIME_STAMP = datetime.datetime.now().strftime('%H-%M-%S')
        self.chip_list = chips
        self.WHOLE_NAME = '_'.join(self.chip_list)
        self.LOG_FILE_PATH = self.WHOLE_NAME + '_' +  self.TIME_STAMP + '.log'
        self.LOG_FILE = open(self.LOG_FILE_PATH, 'w', encoding='utf-8')
        self.SHELL_LOG_FILE = open('SHELL_LOG_' + self.WHOLE_NAME + '_' +  self.TIME_STAMP + '.log', 'w', encoding='utf-8')
        self.OUTPUT_PATH = self.WHOLE_NAME + '_' +  self.TIME_STAMP + '/'
        self.RESULT_FILE = ''
        if os.path.exists(self.OUTPUT_PATH):
            pass
        else:
            os.makedirs(self.OUTPUT_PATH) 

    # 写入一行日志文件, 自动添加换行符
    def write_log(self, content):
        self.LOG_FILE.writelines(content + '\n')
    
    # 读取文件至 ndarray
    def read_file_nd(self, filePath, head=0, raw=False):
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
            line = list(map(str, line.split()))
            data.append(line)
        if raw:
            pass
        else:
            data = np.array(data)
        time.sleep(0.5)
        return data

    # 装饰器, 输出函数运行时间并写入日志
    def time_indicator(func):

        @wraps(func)
        def wrapper(self, *args, **kw):
            print('-------------------------')
            self.write_log('-------------------------')
            start_time = str(datetime.datetime.now())
            print(func.__name__ + ' start at ' + start_time)
            self.write_log(func.__name__ + ' start at ' + start_time)
            res = func(self, *args, **kw)
            end_time = str(datetime.datetime.now())
            print(func.__name__ + ' end at ' + start_time)
            self.write_log(func.__name__ + ' end at ' + end_time)
            return res
        return wrapper

    # 执行 shell 命令并写入日志
    def run_shell(self, command, fileout=True): 
        command = command.replace('$dir', self.OUTPUT_PATH)
        start_time = datetime.datetime.now().strftime('%y-%m-%d %I:%M:%S %f')
        self.write_log('-------------------------')
        self.write_log('Run command > ' + command)
        self.write_log('Start at ' + start_time)
        if fileout:
            subprocess.run(args=command, shell=True, check=True, stdout=self.SHELL_LOG_FILE)
        else:
            subprocess.run(args=command, shell=True, check=True, stdout=self.SHELL_LOG_FILE)
        end_time = datetime.datetime.now().strftime('%y-%m-%d %I:%M:%S %f')
        self.write_log('Finished at ' + end_time)

    # 查找反义链及无效链
    @time_indicator
    def find_reverse(self, file_name):
        f = open(file_name + '.bim', encoding='utf-8')
        lines = f.readlines()
        f.close()
        result_list = []
        for i in range(len(lines)):
            line = list(map(str, lines[i].split()))
            result_list.append([line[1], line[0], line[3], line[4], line[5], line[2]])

        SNP_DICT = {
            'A': 'T',
            'C': 'G',
            'T': 'A',
            'G': 'C',
        }
        reverse_list = [
            ['chromosome', '0_idx_position', 'snp_name', 'genetic_distance', 'allele_1', 'allele_2', 'reference',
            'reference_rev', 'strand']]

        print('Reversing...', file_name + '.bim')
        time.sleep(0.5)
        for i in trange(len(result_list)):
            strand = 'ambiguous'
            ID, CHR, POS, A1, A2, DIS = result_list[i]

            # 跳过染色体号和位置信息不正确的位点, 跳过性染色体, 将其直接置为 ambiguous
            if POS == '0' or CHR == 'X' or CHR == 'Y' or int(CHR) > 18 or int(CHR) < 1:
                reverse_list.append([CHR, POS, ID, DIS, A1, A2, '', '', strand])
                continue
            cmd_samtools = self.SAMTOOLS_PATH + ' faidx ' + self.PIG_FA_PATH + ' chr' + CHR + ':' + str(POS) + '-' + str(POS)
            status, cmd_output = subprocess.getstatusoutput(cmd_samtools)

            if status == 0:
                A1_ = '' + A1
                A2_ = '' + A2
                X = ''.join(cmd_output.split('\n')[1:]).upper()

                # samtools 调用结果不正确的直接置为 ambiguous
                if X not in SNP_DICT:
                    reverse_list.append([CHR, POS, ID, DIS, A1, A2, '', '', strand])
                    continue


                # 两个位点都是0：ambiguous
                if A1_ == '0' and A2_ == '0' and X != 'N':
                    strand = 'ambiguous'
                    reverse_list.append([CHR, POS, ID, DIS, A1, A2, X, SNP_DICT[X], strand])
                elif A1_ == 'N' or A2_ == 'N' or X == 'N':
                    strand = 'forward'
                    reverse_list.append([CHR, POS, ID, DIS, A1, A2, X, X, strand])
                else:
                    # 将带 0 的处理成纯合位点
                    if A1 == '0' and A2 != '0':
                        A1_ = A2
                        A2_ = A2
                    if A2 == '0' and A1 != '0':
                        A1_ = A1
                        A2_ = A1
                    A = A1_ + A2_

                    # 杂合位点
                    if A1_ != A2_:
                        # 有REF，判断为正义
                        if X in A:
                            strand = 'forward'
                        # 无REF，判断为反义
                        else:
                            strand = 'reverse'

                    # 纯合位点
                    else:
                        # 既不是REF也不是REF的互补，判断为正义
                        if X != A1_ and X != SNP_DICT[A1_]:
                            strand = 'forward'
                        # 是REF，判断为正义
                        if X == A1_:
                            strand = 'forward'
                        # 是REF的互补，判断为反义
                        if X == SNP_DICT[A1_]:
                            strand = 'reverse'
                    reverse_list.append([CHR, POS, ID, DIS, A1, A2, X, SNP_DICT[X], strand])
            else:
                reverse_list.append([CHR, POS, ID, DIS, A1, A2, '', '', strand])

        # 结果写入文件
        output = open(self.OUTPUT_PATH + file_name + '_reverse.txt', 'w', encoding='utf-8')
        output1 = open(self.OUTPUT_PATH + file_name + '_reverse_ID.txt', 'w', encoding='utf-8')
        output2 = open(self.OUTPUT_PATH + file_name + '_ambiguous_ID.txt', 'w', encoding='utf-8')
        for reverse in reverse_list:
            output.writelines(' '.join(reverse) + '\n')
            if reverse[-1] == 'reverse':
                output1.writelines(reverse[2] + '\n')
            if reverse[-1] == 'ambiguous':
                output2.writelines(reverse[2] + '\n')
        output.close()
        output1.close()
        output2.close()

    # 转换 reverse 链, 去除 ambiguous 链 
    @time_indicator
    def flip(self, file_name):
        cmd = self.PLINK_PATH + ' --bfile $chip --flip $dir$chip_reverse_ID.txt --exclude $dir$chip_ambiguous_ID.txt --make-bed --out $dir$chip_correct'.replace('$chip', file_name)
        self.run_shell(cmd)

        f = open('$dir$chip_correct.bim'.replace('$chip', file_name).replace('$dir', self.OUTPUT_PATH), encoding='utf-8')
        lines = f.readlines()
        GENE_NUM = len(lines)
        f.close()
        f = open('$dir$chip_correct.fam'.replace('$chip', file_name).replace('$dir', self.OUTPUT_PATH), encoding='utf-8')
        lines = f.readlines()
        SAMPLE_NUM = len(lines)
        f.close()

        self.write_log('Data < ' + file_name + ' > corrected!')
        self.write_log('Data < ' + file_name + ' > has ' + str(GENE_NUM) + ' variants and ' + str(SAMPLE_NUM) + ' individuals before imputation.')

    # 根据数据库的重复位点对应关系进行名称替换
    @time_indicator
    def replace_common_id(self, file_name):
        cmd = self.PLINK_PATH + ' --bfile $dir$chip_correct --update-map duplicate_2.txt --update-name --make-bed --out $dir$chip_correct_new'.replace('$chip', file_name)
        self.run_shell(cmd)

    # 合并所有处理后的芯片文件
    @time_indicator
    def binary_merge(self):
        chip_list = self.chip_list
        for i in range(len(chip_list) - 1):
            chip1 = chip_list[i]
            chip2 = chip_list[i + 1]
            if i == 0:
                cmd1 = self.PLINK_PATH + ' --bfile $dir$chip1_correct_new --bmerge $dir$chip2_correct_new.bed $dir$chip2_correct_new.bim $dir$chip2_correct_new.fam --recode vcf --out $dir$chip1_$chip2_joint'.replace('$chip1', chip1).replace('$chip2', chip2)
                try:
                    self.run_shell(cmd1)
                except:
                    print('Missnp caught')
                    self.write_log('Missnp caught')
                    cmd_error = self.PLINK_PATH + ' --bfile $dir$chip2_correct_new --exclude $dir$chip1_$chip2_joint.missnp --make-bed --out $dir$chip2_correct_new'.replace('$chip1', chip1).replace('$chip2', chip2)
                    self.run_shell(cmd_error)
                    self.run_shell(cmd1)

            else:
                chip0 = chip_list[i - 1]
                cmd2 = self.PLINK_PATH + ' --bfile $dir$chip0_$chip1_joint --bmerge $dir$chip2_correct_new.bed $dir$chip2_correct_new.bim $dir$chip2_correct_new.fam --recode vcf --out $dir$chip1_$chip2_joint'.replace('$chip0', chip0).replace('$chip1', chip1).replace('$chip2', chip2)

                try:
                    self.run_shell(cmd2)
                except:
                    self.write_log('Missnp caught')
                    cmd_error = self.PLINK_PATH + ' --bfile $dir$chip2_correct_new --exclude $dir$chip1_$chip2_joint.missnp --make-bed --out $dir$chip2_correct_new'.replace('$chip1', chip1).replace('$chip2', chip2)
                    self.run_shell(cmd_error)
                    self.run_shell(cmd2)

        chip1 = self.chip_list[-2]
        chip2 = self.chip_list[-1]
        cmd_last = self.PLINK_PATH + ' --bfile $dir$chip1_$chip2_joint --recode vcf --out $dir$whole_name'.replace('$chip1', chip1).replace('$chip2', chip2).replace('$whole_name', self.WHOLE_NAME)
        self.run_shell(cmd_last)

    # beagle 填充
    @time_indicator
    def beagle_chip_list(self):
        chip_list = self.chip_list

        # beagle 命令的输出直接在命令行进行, 方便查看进度
        cmd = 'java -jar ' + self.BEAGLE_PATH + ' gt=$dir$whole_name.vcf out=$dir$whole_name_impute'.replace('$whole_name', self.WHOLE_NAME)
        self.run_shell(cmd, fileout=False)

        gzip_cmd = 'gunzip $dir$whole_name_impute.vcf.gz'.replace('$whole_name', self.WHOLE_NAME)
        self.run_shell(gzip_cmd)

        plink_cmd = self.PLINK_PATH + ' --vcf $dir$whole_name_impute.vcf --recode A --out $dir$whole_name_impute'.replace('$whole_name', self.WHOLE_NAME)
        self.run_shell(plink_cmd)

        recode_f = self.read_file_nd('$dir$whole_name_impute.raw'.replace('$whole_name', self.WHOLE_NAME).replace('$dir', self.OUTPUT_PATH), head=1, raw=True)
        row, col = len(recode_f), len(recode_f[0])
        now = datetime.datetime.now()
        year, month, day, hour, minute = str(now.year), str(now.month), str(now.day), str(now.hour), str(now.minute)
        outputfile_path = 'pi.' + str(col - 6) + '.' + str(row) + '.' + ''.join([year, month, day, hour, minute]) + '.geno_joint.txt'
        self.RESULT_FILE = outputfile_path
        output = open(outputfile_path, 'w', encoding='utf-8')
        for line in recode_f:
            output.writelines(' '.join(line[1:2] + line[6:]) + '\n')
        output.close()
        self.write_log(str(col - 6) + ' variants and ' + str(row) + ' individuals remain after imputation.')

    def run(self):
        chip_list = self.chip_list
        start_time = datetime.datetime.now().strftime('%y-%m-%d %I:%M:%S %f')
        self.write_log('Imputation start at ' + start_time)
        
        # 操作流程
        for cp in chip_list:
            self.find_reverse(cp)
            self.flip(cp)
            self.replace_common_id(cp)
        self.binary_merge()
        self.beagle_chip_list()

        end_time = datetime.datetime.now().strftime('%y-%m-%d %I:%M:%S %f')
        self.write_log('Imputation end at ' + end_time)
        self.LOG_FILE.close()
        self.SHELL_LOG_FILE.close()
        print('Log file', self.LOG_FILE_PATH)
        print('Result file', self.RESULT_FILE)

        # 将结果文件名返回
        return self.RESULT_FILE


if __name__ == '__main__':

    cp_l = sys.argv[1:]
    impute = Union(cp_l)
    impute.run()
