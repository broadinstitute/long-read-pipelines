import sys
import re
import subprocess
import hashlib

"""
Given a GCS path, detect what kind of data we're being given and return a single unmapped BAM file.
"""

bams = []
fqs = []

for arg in sys.argv[1:]:
    if bool(re.search('(\.f(ast)?q)(\.gz)?$', arg)):
        fqs.append(arg)
    elif bool(re.search('\.bam$', arg)):
        bams.append(arg)

if len(fqs) > 0:
    # Fastq files don't have a read group ID to auto-detect, so make one up from a sha1 code of the input files.
    hc = hashlib.sha1(",".join(fqs).encode()).hexdigest()[0:8]
    header = "\t".join(['@HD', 'VN:1.5', 'SO:unknown'])
    rg = "\t".join(['@RG', 'ID:' + hc, 'SM:unknown'])

    # This rather complicated one-liner manually transforms a collection of fastq files into an unmapped BAM file.
    cmd = f'((echo "{header}") && (echo "{rg}") && (zcat {" ".join(fqs)} | paste - - - - | sed "s/^@//" | awk -F "\\t" \'{{ a = sub(/ .+/, "", $1); print $a, "4", "*", "0", "255", "*", "*", "0", "0", $2, $4, "RG:Z:{hc}" }}\' | sed "s/ /\\t/g")) | samtools view -b > unmapped.bam'
elif len(bams) > 1:
    # Merge multiple unmapped BAM files into a single unmapped BAM file.
    cmd = f'samtools cat -o unmapped.bam {" ".join(bams)}'
elif len(bams) == 1:
    # Just move the single BAM file we were given to the name we expect.
    cmd = f'mv {arg} unmapped.bam'

print(cmd)
subprocess.call(cmd, shell=True)
