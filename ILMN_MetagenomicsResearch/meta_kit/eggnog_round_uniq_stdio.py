#!/usr/bin/env python

import sys


out_dict = dict()
for line in sys.stdin.readlines():
    name, number = line.strip().split('\t')
    if name in out_dict.keys():
        out_dict[name] += float(number)
    else:
        out_dict[name] = float(number)

for k,v in out_dict.items():
    # print(k, round(v), sep='\t')
    sys.stdout.write('{}\t{}\n'.format(k,round(v)))
