#!/usr/bin/env python
#_*_coding:utf-8_*_

from collections import Counter
import numpy as np
import re
import math

def Count_1(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum
def Count_2(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 0 for i in cutoffNums]
    code = []
    for cutoff in cutoffNums:
        myCount = 0
        if cutoff == 0:
            code.append(0)
        else:
            for i in range(len(sequence)):
                if sequence[i] in aaSet:
                    myCount += 1
                    if myCount == cutoff:
                        code.append((i + 1) / len(sequence) * 100)
                        break
            if myCount == 0:
                code.append(0)
    return code

AA = 'ACDEFGHIKLMNPQRSTVWY'
kw = {'order': 'ACDEFGHIKLMNPQRSTVWY'}

def get_features(fasta):
    def APAAC(lambdaValue=2, w=0.05, **kw):
        with open("./File_list/PAAC.txt", "r") as f:
            records = f.readlines()
        AA = ''.join(records[0].rstrip().split()[1:])
        AADict = {}
        for i in range(len(AA)):
            AADict[AA[i]] = i
        AAProperty = []
        AAPropertyNames = []
        for i in range(1, len(records) - 1):
            array = records[i].rstrip().split() if records[i].rstrip() != '' else None
            AAProperty.append([float(j) for j in array[1:]])
            AAPropertyNames.append(array[0])
            AAProperty1 = []
        for i in AAProperty:
            meanI = sum(i) / 20
            fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
            AAProperty1.append([(j - meanI) / fenmu for j in i])
        encodings = []

        for i in fasta:
            nam ,sequence = i[0],re.sub('-', '', i[1])
            code = []
            theta = []
            for n in range(1, lambdaValue + 1):
                for j in range(len(AAProperty1)):
                    theta.append(sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in range(len(sequence) - n)]) / (len(sequence) - n))
            myDict = {}
            for aa in AA:
                myDict[aa] = sequence.count(aa)
            code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
            code = code + [w * value / (1 + w * sum(theta)) for value in theta]
            encodings.append(code)
        return encodings
    def DPC(**kw):
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        encodings = []


        AADict = {}
        for i in range(len(AA)):
            AADict[AA[i]] = i

        for i in fasta:
            name,sequence=i[0], re.sub('-', '', i[1])
            code = []
            tmpCode = [0] * 400
            for j in range(len(sequence) - 2 + 1):
                tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
            if sum(tmpCode) != 0:
                tmpCode = [i/sum(tmpCode) for i in tmpCode]
            code = code + tmpCode
            encodings.append(code)
        return encodings

    def Rvalue(aa1, aa2, AADict,Matrix):
        return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

    def PAAC():
        dataFile = r'./File_list/PAAC.txt'
        with open(dataFile) as f:
            record = f.readlines()
        AA = ''.join(record[0].rstrip().split()[1:])  # 读取氨基酸序列可知AA = "ARNDCQEGHILKMFPSTWYV"
        AADict = {}  # 用来存放氨基酸在AA中的位置
        for i in range(len(AA)):
            AADict[AA[i]] = i
        AAproperty = []  # 用来提取PAAC相关属性
        AApropertyName = []
        for i in range(1, len(record)):  # 从第二行开始，第一行为氨基酸序列，第二行开始为属性
            array = record[i].rstrip().split() if record[i].rstrip().split() != '' else None
            AAproperty.append([float(j) for j in array[1:]])
            AApropertyName.append(array[0])
        AAproperty1 = []
        for i in AAproperty:
            meanI = sum(i) / 20  # 每一行即每个属性对应氨基酸的均值
            fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)  # 每一行的标准差
            AAproperty1.append([(j - meanI) / fenmu for j in i])  # z-score
        encoding = []
        for i in fasta:
            name, sequences = i[0], re.sub('-', '', i[1])
            code = []
            theta = []
            for n in range(1, 3):
                theta.append(sum([Rvalue(sequences[j], sequences[j + n], AADict, AAproperty1) for j in
                                  range(len(sequences) - n)]) / (len(sequences) - n))
            myDict = {}
            for aa in AA:
                myDict[aa] = sequences.count(aa)
            code = code + [myDict[aa] / (1 + 0.05 * sum(theta)) for aa in AA]  # 计入特征， 20 维
            code = code + [(0.05 * j) / (1 + 0.05 * sum(theta)) for j in theta]  # 计入特征 2 维 故共22 维特征
            encoding.append(code)
        return encoding


    def ASDC():
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)
        for i in fasta:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            sum = 0
            pair_dict = {}
            for pair in aaPairs:
                pair_dict[pair] = 0
            for j in range(len(sequence)):
                for k in range(j + 1, len(sequence)):
                    if sequence[j] in AA and sequence[k] in AA:
                        pair_dict[sequence[j] + sequence[k]] += 1
                        sum += 1
            for pair in aaPairs:
                code.append(pair_dict[pair] / sum)
            encodings.append(code)
        return encodings

    def QSOrder(nlag=2, w=0.05, **kw):
        AA = 'ACDEFGHIKLMNPQRSTVWY'
        AA1 = 'ARNDCQEGHILKMFPSTWYV'

        DictAA = {}
        for i in range(len(AA)):
            DictAA[AA[i]] = i

        DictAA1 = {}
        for i in range(len(AA1)):
            DictAA1[AA1[i]] = i
        with open("./File_list/Schneider-Wrede.txt", "r") as f:
            records = f.readlines()[1:]
        AADistance = []
        for i in records:
            array = i.rstrip().split()[1:] if i.rstrip() != '' else None
            AADistance.append(array)
        AADistance = np.array(
            [float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape(
            (20, 20))

        with open("./File_list/Grantham.txt", "r") as f:
            records = f.readlines()[1:]
        AADistance1 = []
        for i in records:
            array = i.rstrip().split()[1:] if i.rstrip() != '' else None
            AADistance1.append(array)
        AADistance1 = np.array(
            [float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
            (20, 20))

        encodings = []

        for i in fasta:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = []
            arraySW = []
            arrayGM = []
            for n in range(1, nlag + 1):
                arraySW.append(
                    sum([AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in
                         range(len(sequence) - n)]))
                arrayGM.append(sum(
                    [AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in
                     range(len(sequence) - n)]))
            myDict = {}
            for aa in AA1:
                myDict[aa] = sequence.count(aa)
            for aa in AA1:
                code.append(myDict[aa] / (1 + w * sum(arraySW)))
            for aa in AA1:
                code.append(myDict[aa] / (1 + w * sum(arrayGM)))
            for num in arraySW:
                code.append((w * num) / (1 + w * sum(arraySW)))
            for num in arrayGM:
                code.append((w * num) / (1 + w * sum(arrayGM)))
            encodings.append(code)
        return encodings
    print('Feature extraction...')
    encoding = []
    encoding.append(APAAC())
    encoding.append(DPC())
    encoding.append(ASDC())
    encoding.append(PAAC())
    encoding.append(QSOrder())
    return np.column_stack(encoding)

