#!/usr/bin/env bash
# Reproduce the bcf2hmmibd equivalence proof.
#
#   verify_equivalence.sh [bindir]
#
# bindir (default /usr/local/bin) must contain `hmmibd-rs` and `bcf2hmmibd`;
# `hmmIBD` (the bundled C tool) is used for the optional second comparison if
# present. `bcftools` and `python3` must be on PATH.
#
# Proves:
#   1) bcf2hmmibd table + `hmmibd-rs -i`  ==  `hmmibd-rs --from-bcf`  (byte-identical)
#   2) C hmmIBD  ==  hmmibd-rs  on the same table                      (numerically identical)
set -euo pipefail

BINDIR="${1:-/usr/local/bin}"
HMMIBD_RS="$BINDIR/hmmibd-rs"
BCF2="$BINDIR/bcf2hmmibd"
CHMM="$BINDIR/hmmIBD"
HERE="$(cd "$(dirname "$0")" && pwd)"

command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH"; exit 2; }
[ -x "$HMMIBD_RS" ] || { echo "ERROR: $HMMIBD_RS not found/executable"; exit 2; }
[ -x "$BCF2" ]      || { echo "ERROR: $BCF2 not found/executable"; exit 2; }

WORK="$(mktemp -d)"
cd "$WORK"

python3 "$HERE/gen_test_vcf.py"
bcftools view -Ob -o test.bcf test.vcf
bcftools index test.bcf

"$HMMIBD_RS" -i test.bcf --from-bcf -o A >/dev/null 2>&1
"$BCF2"      -i test.bcf -o table.txt      2>/dev/null
"$HMMIBD_RS" -i table.txt -o C          >/dev/null 2>&1

echo
echo "== (1) bcf2hmmibd table + 'hmmibd-rs -i'  vs  'hmmibd-rs --from-bcf' =="
if diff -q A.hmm.txt C.hmm.txt >/dev/null && diff -q A.hmm_fract.txt C.hmm_fract.txt >/dev/null; then
    echo "PASS: byte-identical (md5 $(md5sum < A.hmm.txt | cut -c1-12).. / $(md5sum < A.hmm_fract.txt | cut -c1-12)..)"
else
    echo "FAIL: outputs differ"; diff A.hmm.txt C.hmm.txt | head; exit 1
fi

if [ -x "$CHMM" ]; then
    echo
    echo "== (2) C hmmIBD  vs  hmmibd-rs  on the same bcf2hmmibd table =="
    "$CHMM" -i table.txt -o Cc >/dev/null 2>&1
    python3 "$HERE/cmp.py" Cc.hmm_fract.txt C.hmm_fract.txt
    python3 "$HERE/cmp.py" Cc.hmm.txt       C.hmm.txt
    echo "(only cosmetic log_p exponent formatting differs: C '-3.8e+01' vs Rust '-3.8e1')"
fi

echo
echo "workdir: $WORK"
