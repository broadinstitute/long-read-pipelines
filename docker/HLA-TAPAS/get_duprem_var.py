#!/usr/bin/env python
import sys
INPUT = sys.argv[1]

RS_missing = {}

with open(INPUT + ".lmiss") as f:
	for line in f:
		if line.find("CHR") == -1:
			line = line.rstrip().split()
			RS_missing[line[1]] = float(line[4])


OUT = open(INPUT + ".remdup.snp","w")

pos_snp_dic = {}

with open(INPUT + ".bim") as f:
	for line in f:
		line = line.rstrip().split()
		pos_snp_dic.setdefault((line[0],line[3],line[4],line[5]),[]).append(line[1])

for pos in pos_snp_dic:
	if len(pos_snp_dic[pos]) > 1:
		snp_list = pos_snp_dic[pos]
		missings = {}
		for snp in snp_list:
			missings[snp] = RS_missing[snp]
		name, min_mis = min(missings.items(), key = lambda x: x[1])
		for snp in snp_list:
			if snp != name:
				print(snp, file = OUT)

OUT.close()
