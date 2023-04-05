#!/usr/bin/env python
# _*_coding:utf-8_*_

import re

def check_sequences(data,fastas_seq):
    seq_len = []
    for i in fastas_seq:
        seq = re.sub('-', '', i[1])
        seq_len.append(len(seq.rstrip()))
        print(seq_len)
    if re.search('>', data) == None:
        return False
    elif min(seq_len) == 0:
        return False
    else:
        return True










