import gzip
import sys

MINIMUM_STR_LENGTH = 3
VALID_NUCLEOTIDES = {'A', 'C', 'G', 'T'}

def compute_repeat_unit(nuc_sequence):
    """Takes a nucleotide sequence and checks whether it consists of multiple exact repeats of a smaller sub-sequence.
    If yes, it returns this sub-sequence, otherwise it returns None.

    Examples:
        ACACT returns None
        ACACA returns AC
        ACTACTA returns ACT
    """
    for repeat_size in range(1, len(nuc_sequence)//2 + 1 + 1):
        repeat_unit = nuc_sequence[:repeat_size]
        repeats = repeat_unit * (1 + len(nuc_sequence)//repeat_size)
        if nuc_sequence == repeats[:len(nuc_sequence)]:
            return repeat_unit

    return None


#compute_repeat_unit("ACAC")

input_vcf = sys.argv[1]

fopen = gzip.open if input_vcf.endswith(".gz") else open

for line in fopen(input_vcf):
    if line.startswith("#"):
        sys.stdout.write(line)
        continue
    
    fields = line.split("\t")
    ref = fields[3].upper()
    alt = fields[4].upper()    
    #repeat_unit = compute_repeat_unit(alt if len(alt) >= MINIMUM_STR_LENGTH else ref)
    repeat_unit = compute_repeat_unit(alt[1:] if len(alt[1:]) >= MINIMUM_STR_LENGTH else ref[1:])
    if repeat_unit is None:
        continue

    info_field = fields[7]
    if "RU" not in info_field and repeat_unit is not None:
        fields[7] = "RU={};{}".format(repeat_unit, info_field)

    id_field = ""
    if repeat_unit:
        id_field += "RU{}:{}:".format(len(repeat_unit), repeat_unit)

    if len(ref) != len(alt):
        if len(ref) > len(alt):
            id_field += "DEL:"
        elif len(ref) < len(alt):
            id_field += "INS:"

        id_field += "{}=>{}".format(len(ref)-1, len(alt)-1)

    fields[2] = id_field

    line = "\t".join(fields)

    sys.stdout.write(line)
