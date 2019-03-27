import sys
import gzip
import collections

rhash = collections.OrderedDict()
total = 0

for report in sys.argv[1:]:
    fopen = gzip.open if report.endswith(".gz") else open

    for line in fopen(report):
        if not line.startswith("ZMW") and "," in line:
            fields = line.split(",")
            key = fields[0]
            count = int(fields[1])

            if not key in rhash:
                rhash[key] = 0

            rhash[key] += count
            total += count

sys.stdout.write("ZMW Yield\n")
for key in rhash:
    sys.stdout.write("%s,%d,%.2f%%\n" % (key, rhash[key], 100.0*rhash[key]/total))

sys.stdout.write("\n\n")
