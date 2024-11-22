#!/use/bin/env python
import sys
REFBIMFILE = sys.argv[1]
OUTBIMFILE = sys.argv[2]

REF = {}
with open(REFBIMFILE) as f:
	for line in f:
		line = line.rstrip().split()
		CHRPOS = ":".join([line[0],line[3]])
		VARNAME = line[1]
		A1 = line[4]
		A2 = line[5]
		REF.setdefault(CHRPOS,[]).append([VARNAME,A1,A2])
OUT = open(OUTBIMFILE, "w")
with open(OUTBIMFILE + ".old") as f:
	for line in f:
		line = line.rstrip().split()
		CHRPOS = ":".join([line[0],line[3]])
		VARNAME = line[1]
		A1 = line[4]
		A2 = line[5]
		if CHRPOS in REF:
			for candidate in REF[CHRPOS]:
				if (A1 == candidate[1] and A2 == candidate[2]) or (A1 == candidate[2] and A2 == candidate[1]):
					VARNAME = candidate[0]
		out = " ".join([line[0], VARNAME, line[2], line[3], line[4], line[5]])
		print(out, file = OUT)
OUT.close()
