#!/usr/bin/env python
# _*_coding:utf-8_*_

import re

def read_protein_sequences(data):
    records = data.split('>')[1:]
    fasta_sequences = []
    sequence_name = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = array[0]
        fasta_sequences.append(([name, sequence]))
        sequence_name.append(str(name))
    return fasta_sequences, sequence_name





