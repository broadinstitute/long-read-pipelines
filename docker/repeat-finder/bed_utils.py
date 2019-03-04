from collections import namedtuple
import gzip
import os
import re


BedRecord = namedtuple('BedRecord', ['chrom', 'start', 'end', 'repeat_unit', 'repeat_count'])


def parse_bed_file(bed_path, limit=None):

    fopen = gzip.open if bed_path.endswith("gz") else open

    with fopen(os.path.expanduser(bed_path), mode="rt") as f:
        for i, line in enumerate(f):

            if limit is not None and i >= limit:
                break

            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            try:
                chrom = re.sub("^chr", "", fields[0]).upper()
                if '_' in chrom or len(chrom) > 5:
                    continue

                start_0based = int(fields[1])
                end_1based = int(fields[2])
                repeat_unit = None
                repeat_count = None

                if len(fields) > 3 and len(set(fields[3]) - {'A', 'C', 'G', 'T', 'N'}) == 0:
                    repeat_unit = fields[3].upper()
                    repeat_count = float(fields[4])
                elif len(fields) > 4 and len(set(fields[4]) - {'A', 'C', 'G', 'T', 'N'}) == 0:
                    repeat_count = float(fields[3])
                    repeat_unit = fields[4].upper()
                elif len(fields) > 3:
                    raise ValueError(f"column[4] is not a nucleotide sequence: {repeat_unit}")

            except Exception as e:
                raise ValueError(f"Unable to parse line {i}: {e}\n{line}")

            yield BedRecord(chrom, start_0based, end_1based, repeat_unit, repeat_count)
