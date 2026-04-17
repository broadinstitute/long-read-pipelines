#!/usr/bin/env python3

import random
import subprocess
import tempfile
from pathlib import Path


REFERENCE = Path("/Users/shuang/Ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz")
FAI = Path(str(REFERENCE) + ".fai")
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_BAM = REPO_ROOT / "tests" / "fixtures" / "test.dev.ONT.withHardClips.usingGRCh38seq.bam"
SEED = 20260416

IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(seq: str) -> str:
    return seq.translate(IUPAC_COMPLEMENT)[::-1]


def read_fai() -> list[tuple[str, int]]:
    contigs = []
    with FAI.open() as fh:
        for line in fh:
            name, length, *_ = line.rstrip().split("\t")
            if "_" in name:
                continue
            contigs.append((name, int(length)))
    return contigs


def fetch_region_sequence(contig: str, start_1based: int, end_1based: int) -> str:
    region = f"{contig}:{start_1based}-{end_1based}"
    result = subprocess.run(
        ["samtools", "faidx", str(REFERENCE), region],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = result.stdout.splitlines()
    return "".join(lines[1:]).upper()


def mutate_sequence(seq: str, rng: random.Random, rate: float = 0.012) -> str:
    bases = "ACGT"
    seq_list = list(seq)
    mutations = max(1, int(len(seq_list) * rate))

    for idx in rng.sample(range(len(seq_list)), mutations):
        ref = seq_list[idx]
        if ref in bases:
            choices = [b for b in bases if b != ref]
            seq_list[idx] = rng.choice(choices)
        elif ref == "N":
            seq_list[idx] = rng.choice(bases)

    return "".join(seq_list)


def random_ont_quality(length: int, rng: random.Random) -> str:
    return "".join(chr(rng.randint(8, 19) + 33) for _ in range(length))


def build_primary_cigar(read_len: int, rng: random.Random) -> str:
    insertion_len = rng.randint(10, 40)
    deletion_len = rng.randint(6, 28)
    insertion_len_2 = rng.randint(6, 26)
    m_total = read_len - insertion_len - insertion_len_2
    m1 = max(700, int(m_total * rng.uniform(0.20, 0.30)))
    m2 = max(700, int(m_total * rng.uniform(0.24, 0.36)))
    m3 = m_total - m1 - m2

    if m3 < 700:
        deficit = 700 - m3
        m3 = 700
        m2 -= deficit

    return f"{m1}M{insertion_len}I{m2}M{deletion_len}D{m3}M{insertion_len_2}I"


def build_hardclipped_cigar(full_len: int, start: int, seg_len: int, rng: random.Random) -> str:
    left_h = start
    right_h = full_len - start - seg_len
    insertion_len = rng.randint(5, 20)
    deletion_len = rng.randint(3, 16)
    m_total = seg_len - insertion_len
    m1 = max(300, int(m_total * rng.uniform(0.38, 0.58)))
    m2 = m_total - m1

    if m2 < 300:
        deficit = 300 - m2
        m2 = 300
        m1 -= deficit

    ops = []
    if left_h:
        ops.append(f"{left_h}H")
    ops.append(f"{m1}M")
    ops.append(f"{insertion_len}I")
    ops.append(f"{deletion_len}D")
    ops.append(f"{m2}M")
    if right_h:
        ops.append(f"{right_h}H")
    return "".join(ops)


def orient(seq: str, qual: str, is_reverse: bool) -> tuple[str, str]:
    if not is_reverse:
        return seq, qual
    return reverse_complement(seq), qual[::-1]


def make_record(
    qname: str,
    flag: int,
    rname: str,
    pos: int,
    mapq: int,
    cigar: str,
    seq: str,
    qual: str,
    tags: list[str],
) -> str:
    fields = [
        qname,
        str(flag),
        rname,
        str(pos),
        str(mapq),
        cigar,
        "*",
        "0",
        "0",
        seq,
        qual,
    ]
    fields.extend(tags)
    return "\t".join(fields)


def pick_valid_window(
    contig: str, contig_len: int, read_len: int, seed_offset: int, rng: random.Random
) -> tuple[int, str]:
    max_start = contig_len - read_len - 1
    candidate_starts = [
        100_000 + seed_offset,
        max(1, int(contig_len * 0.12) + seed_offset),
        max(1, int(contig_len * 0.36) + seed_offset),
        max(1, int(contig_len * 0.62) + seed_offset),
    ]

    for start in candidate_starts:
        start = min(start, max_start)
        seq = fetch_region_sequence(contig, start, start + read_len - 1)
        if len(seq) == read_len and seq.count("N") / read_len < 0.10:
            return start, seq

    for _ in range(30):
        start = rng.randint(1, max_start)
        seq = fetch_region_sequence(contig, start, start + read_len - 1)
        if len(seq) == read_len and seq.count("N") / read_len < 0.10:
            return start, seq

    raise RuntimeError(f"unable to fetch a suitable window for {contig}:{read_len}")


def main() -> None:
    if not REFERENCE.exists():
        raise FileNotFoundError(REFERENCE)
    if not FAI.exists():
        raise FileNotFoundError(FAI)

    rng = random.Random(SEED)
    read_lengths = [12940, 13680, 14220, 14780, 15110, 15360, 15740, 16190, 16850, 18130]
    contigs = read_fai()
    canonical = contigs[:25]

    if OUT_BAM.exists():
        OUT_BAM.unlink()

    header_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    header_lines.extend(f"@SQ\tSN:{name}\tLN:{length}" for name, length in canonical)
    header_lines.append("@PG\tID:generate_test_bam.usingGRCh38seq.py\tPN:generate_test_bam.usingGRCh38seq.py\tVN:1.0")

    records = []

    for i, read_len in enumerate(read_lengths, start=1):
        primary_ref, primary_ref_len = canonical[(i * 2) % 24]
        primary_start, ref_seq = pick_valid_window(primary_ref, primary_ref_len, read_len, i * 5000, rng)
        full_seq = mutate_sequence(ref_seq, rng)
        full_qual = random_ont_quality(read_len, rng)

        qname = f"simGRCh38ONT_{i:04d}"
        primary_reverse = (i % 3 == 0)
        primary_flag = 16 if primary_reverse else 0
        primary_seq, primary_qual = orient(full_seq, full_qual, primary_reverse)
        primary_cigar = build_primary_cigar(read_len, rng)

        sa_ref, sa_ref_len = canonical[(i * 2 + 5) % 24]
        sa_pos = min(1_000_000 + i * 175_000, sa_ref_len - 20_000)
        sa_start = int(read_len * (0.08 + 0.02 * (i % 4)))
        sa_len = int(read_len * (0.34 + 0.03 * (i % 3)))
        sa_len = min(sa_len, read_len - sa_start - 1200)
        sa_len = max(sa_len, 4200)
        sa_segment_seq = full_seq[sa_start : sa_start + sa_len]
        sa_segment_qual = full_qual[sa_start : sa_start + sa_len]
        sa_reverse = (i % 2 == 0)
        sa_flag = 2048 | (16 if sa_reverse else 0)
        sa_seq, sa_qual = orient(sa_segment_seq, sa_segment_qual, sa_reverse)
        sa_cigar = build_hardclipped_cigar(read_len, sa_start, sa_len, rng)

        sec_ref, sec_ref_len = canonical[(i * 2 + 11) % 24]
        sec_pos = min(2_000_000 + i * 135_000, sec_ref_len - 20_000)
        sec_start = int(read_len * (0.22 + 0.015 * (i % 5)))
        sec_len = int(read_len * (0.42 + 0.02 * ((i + 1) % 4)))
        sec_len = min(sec_len, read_len - sec_start - 800)
        sec_len = max(sec_len, 5000)
        sec_segment_seq = full_seq[sec_start : sec_start + sec_len]
        sec_segment_qual = full_qual[sec_start : sec_start + sec_len]
        sec_reverse = (i % 4 in (1, 2))
        sec_flag = 256 | (16 if sec_reverse else 0)
        sec_seq, sec_qual = orient(sec_segment_seq, sec_segment_qual, sec_reverse)
        sec_cigar = build_hardclipped_cigar(read_len, sec_start, sec_len, rng)

        tags = [f"RG:Z:SIM{i % 2 + 1}"]
        records.append(
            (
                primary_ref,
                primary_start,
                make_record(
                    qname,
                    primary_flag,
                    primary_ref,
                    primary_start,
                    60,
                    primary_cigar,
                    primary_seq,
                    primary_qual,
                    tags + [f"tp:A:P", f"qs:i:0", f"qe:i:{read_len}"],
                ),
            )
        )
        records.append(
            (
                sa_ref,
                sa_pos,
                make_record(
                    qname,
                    sa_flag,
                    sa_ref,
                    sa_pos,
                    35,
                    sa_cigar,
                    sa_seq,
                    sa_qual,
                    tags + [f"tp:A:S", f"qs:i:{sa_start}", f"qe:i:{sa_start + sa_len}"],
                ),
            )
        )
        records.append(
            (
                sec_ref,
                sec_pos,
                make_record(
                    qname,
                    sec_flag,
                    sec_ref,
                    sec_pos,
                    22,
                    sec_cigar,
                    sec_seq,
                    sec_qual,
                    tags + [f"tp:A:S", f"qs:i:{sec_start}", f"qe:i:{sec_start + sec_len}"],
                ),
            )
        )

        assert len(full_seq) == read_len
        assert len(primary_seq) == read_len
        assert len(primary_qual) == read_len
        assert len(sa_seq) == sa_len
        assert len(sa_qual) == sa_len
        assert len(sec_seq) == sec_len
        assert len(sec_qual) == sec_len

    records.sort(key=lambda item: (item[0], item[1], item[2].split("\t", 1)[0]))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        sam_path = tmpdir_path / "synthetic.sam"
        unsorted_bam = tmpdir_path / "synthetic.unsorted.bam"

        with sam_path.open("w", encoding="ascii") as sam_file:
            for line in header_lines:
                sam_file.write(line)
                sam_file.write("\n")
            for _, _, record in records:
                sam_file.write(record)
                sam_file.write("\n")

        subprocess.run(
            ["samtools", "view", "-b", "-o", str(unsorted_bam), str(sam_path)],
            check=True,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(OUT_BAM), str(unsorted_bam)],
            check=True,
        )
        subprocess.run(["samtools", "quickcheck", str(OUT_BAM)], check=True)


if __name__ == "__main__":
    main()
