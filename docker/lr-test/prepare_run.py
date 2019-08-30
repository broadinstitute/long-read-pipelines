import sys
import re
import subprocess

bams = []
fqs = []

for arg in sys.argv[1:]:
    if bool(re.search('(\.f(ast)?q)(\.gz)?$', arg)):
        fqs.append(arg)
    elif bool(re.search('\.bam$', arg)):
        bams.append(arg)

if len(fqs) > 0:
    header = "\t".join(['@HD', 'VN:1.5', 'SO:unknown'])
    cmd = f'((echo "{header}") && (zcat {" ".join(fqs)} | paste - - - - | sed "s/^@//" | awk \'{{ print $1, "4", "*", "0", "255", "*", "*", "0", "0", $2, $4 }}\' | sed "s/ /\\t/g")) | samtools view -b > unmapped.bam'
elif len(bams) > 1:
    cmd = f'samtools cat -o unmapped.bam {" ".join(bams)}'
elif len(bams) == 1:
    cmd = f'mv {arg} unmapped.bam'

print(cmd)
subprocess.call(cmd, shell=True)
