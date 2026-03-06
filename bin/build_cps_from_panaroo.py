#!/usr/bin/env python3
"""
Build CPS-focused outputs from Panaroo outputs + serotype reference GenBank.

Approach (sequence-driven, not gene-name-driven):
1) Extract CDS from reference GenBank, excluding pseudogenes/transposon-like entries.
2) BLAST reference CDS against Panaroo pan_genome_reference.fa to map CPS genes to Panaroo clusters.
3) Reconstruct per-isolate gene sequences for mapped clusters, preferring gene_data.csv DNA sequences.
4) For each mapped CPS gene, create an MSA (mafft) across isolates.
5) Concatenate per-gene MSAs in reference order to create per-isolate CPS FASTA and core alignment.
6) Create SNP VCF using snp-sites from concatenated alignment.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

EXCLUDE_TERMS = {
    "transposase",
    "transposon",
    "insertion sequence",
    "integrase",
    "recombinase",
}

PANAROO_METADATA_COLUMNS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}


def run_cmd(cmd: List[str], *, cwd: Path | None = None) -> None:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )


def which_or_raise(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        raise RuntimeError(f"Required executable not found in PATH: {tool}")
    return path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panaroo-dir", required=True, help="Path to Panaroo output directory")
    p.add_argument("--reference-gb", required=True, help="Reference serotype GenBank file")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--min-pident", type=float, default=90.0, help="Minimum BLAST percent identity")
    p.add_argument("--min-qcov", type=float, default=0.8, help="Minimum query coverage (0-1)")
    p.add_argument("--threads", type=int, default=4, help="Threads for BLAST / MAFFT")
    return p.parse_args()


def is_excluded(feature: SeqFeature) -> bool:
    qualifiers = feature.qualifiers
    if "pseudo" in qualifiers:
        return True

    text_fields = []
    for key in ("product", "gene", "note"):
        text_fields.extend(qualifiers.get(key, []))
    text = " ".join(text_fields).lower()
    return any(term in text for term in EXCLUDE_TERMS)


def feature_name(feature: SeqFeature, idx: int) -> str:
    for key in ("locus_tag", "gene", "protein_id"):
        vals = feature.qualifiers.get(key)
        if vals:
            return vals[0]
    return f"cds_{idx:04d}"


def extract_reference_cds(reference_gb: Path, out_fasta: Path) -> List[str]:
    record = SeqIO.read(reference_gb, "genbank")
    kept_names: List[str] = []

    with out_fasta.open("w") as handle:
        i = 0
        for feat in record.features:
            if feat.type != "CDS":
                continue
            i += 1
            if is_excluded(feat):
                continue
            seq = feat.extract(record.seq)
            name = feature_name(feat, i)
            kept_names.append(name)
            handle.write(f">{name}\n{str(seq)}\n")

    if not kept_names:
        raise ValueError("No CPS CDS retained from reference GenBank after filtering")
    return kept_names


def blast_map(ref_cds_fa: Path, pan_ref_fa: Path, threads: int, min_pident: float, min_qcov: float) -> Dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="cps_blast_") as tmpd:
        tmp = Path(tmpd)
        db = tmp / "panaroo_db"
        run_cmd(["makeblastdb", "-in", str(pan_ref_fa), "-dbtype", "nucl", "-out", str(db)])

        out = tmp / "hits.tsv"
        run_cmd(
            [
                "blastn",
                "-query",
                str(ref_cds_fa),
                "-db",
                str(db),
                "-num_threads",
                str(threads),
                "-outfmt",
                "6 qseqid sseqid pident length qlen bitscore",
                "-out",
                str(out),
            ]
        )

        best: Dict[str, Tuple[float, float, str]] = {}
        with out.open() as f:
            for line in f:
                qseqid, sseqid, pident, length, qlen, bitscore = line.strip().split("\t")
                pident_f = float(pident)
                qcov = float(length) / float(qlen)
                bits = float(bitscore)
                if pident_f < min_pident or qcov < min_qcov:
                    continue
                prev = best.get(qseqid)
                if prev is None or bits > prev[0]:
                    best[qseqid] = (bits, qcov, sseqid)

        return {q: v[2] for q, v in best.items()}


def parse_presence_cell(cell: str) -> List[str]:
    if not cell:
        return []
    parts = []
    for token in cell.replace('"', '').split(";"):
        token = token.strip()
        if token:
            parts.append(token)
    return parts


def read_combined_dna_cds(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for rec in SeqIO.parse(path, "fasta"):
        rid = rec.id.strip()
        seq = str(rec.seq)
        seqs[rid] = seq
        desc0 = rec.description.split()[0].strip()
        if desc0 and desc0 not in seqs:
            seqs[desc0] = seq
    return seqs


def read_gene_presence_absence(path: Path) -> Tuple[List[str], Dict[str, Dict[str, List[str]]]]:
    with path.open(newline="") as h:
        reader = csv.DictReader(h)
        fieldnames = reader.fieldnames or []
        isolate_columns = [c for c in fieldnames if c not in PANAROO_METADATA_COLUMNS]

        rows: Dict[str, Dict[str, List[str]]] = {}
        for row in reader:
            gene = (row.get("Gene") or "").strip()
            if not gene:
                continue
            rows[gene] = {
                iso: parse_presence_cell((row.get(iso) or "").strip()) for iso in isolate_columns
            }

    return isolate_columns, rows


def _norm(s: str) -> str:
    return s.strip().lower().replace(" ", "_")


def find_column(fieldnames: List[str], candidates: List[str]) -> Optional[str]:
    if not fieldnames:
        return None
    norm_map = {_norm(f): f for f in fieldnames}
    for c in candidates:
        if _norm(c) in norm_map:
            return norm_map[_norm(c)]
    # soft contains fallback
    for f in fieldnames:
        nf = _norm(f)
        if any(_norm(c) in nf for c in candidates):
            return f
    return None


def read_gene_data(path: Path) -> Tuple[Set[str], Dict[str, Dict[str, List[str]]]]:
    """
    Return: isolates set, mapping panaroo_gene -> isolate -> list[DNA seq].
    """
    with path.open(newline="") as h:
        reader = csv.DictReader(h)
        fieldnames = reader.fieldnames or []

        gene_col = find_column(fieldnames, ["Gene", "gene", "cluster", "gene_id", "panaroo_gene"])
        iso_col = find_column(fieldnames, ["isolate", "genome", "sample", "strain"])
        dna_col = find_column(fieldnames, ["dna_sequence", "dna", "nucleotide_sequence", "sequence", "seq"])

        if not gene_col or not dna_col:
            raise RuntimeError(
                "Could not parse gene_data.csv columns. Need gene and DNA sequence columns. "
                f"Found: {fieldnames}"
            )

        mapping: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
        isolates: Set[str] = set()

        for row in reader:
            gene = (row.get(gene_col) or "").strip()
            dna = (row.get(dna_col) or "").strip()
            if not gene or not dna:
                continue

            isolate = (row.get(iso_col) or "").strip() if iso_col else ""
            if not isolate:
                isolate = "unknown"
            isolates.add(isolate)
            mapping[gene][isolate].append(dna)

    return isolates, mapping


def sanitize_isolate_name(name: str) -> str:
    return name.replace("/", "_").replace(" ", "_")


def write_fasta(path: Path, seqs: Dict[str, str]) -> None:
    with path.open("w") as h:
        for name, seq in seqs.items():
            h.write(f">{name}\n{seq}\n")


def align_gene(in_fa: Path, out_fa: Path, threads: int) -> None:
    which_or_raise("mafft")
    cmd = ["mafft", "--thread", str(threads), "--auto", str(in_fa)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed for {in_fa.name}:\n{proc.stderr}")
    out_fa.write_text(proc.stdout)


def read_alignment(path: Path) -> OrderedDict[str, str]:
    d: OrderedDict[str, str] = OrderedDict()
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = str(rec.seq)
    if not d:
        raise RuntimeError(f"Empty alignment: {path}")
    return d


def best_seq_for_ids(ids: Iterable[str], seq_index: Dict[str, str]) -> str | None:
    best = None
    for gid in ids:
        seq = seq_index.get(gid) or seq_index.get(gid.split("|")[0])
        if seq is None:
            continue
        if best is None or len(seq) > len(best):
            best = seq
    return best


def best_seq_from_list(seqs: Iterable[str]) -> str | None:
    best = None
    for s in seqs:
        if not s:
            continue
        if best is None or len(s) > len(best):
            best = s
    return best


def main() -> int:
    args = parse_args()

    panaroo_dir = Path(args.panaroo_dir)
    output = Path(args.output)
    ref_gb = Path(args.reference_gb)

    pan_ref = panaroo_dir / "pan_genome_reference.fa"
    gpa_csv = panaroo_dir / "gene_presence_absence.csv"
    combined_cds = panaroo_dir / "combined_DNA_CDS.fasta"
    gene_data_csv = panaroo_dir / "gene_data.csv"

    if not pan_ref.exists():
        raise FileNotFoundError(f"Missing required Panaroo file: {pan_ref}")

    which_or_raise("makeblastdb")
    which_or_raise("blastn")
    which_or_raise("snp-sites")

    output.mkdir(parents=True, exist_ok=True)
    cps_gene_align_dir = output / "cps_gene_alignments"
    cps_gene_align_dir.mkdir(exist_ok=True)
    raw_gene_dir = output / "cps_gene_sequences_raw"
    raw_gene_dir.mkdir(exist_ok=True)
    cps_fastas_dir = output / "cps_fastas"
    cps_fastas_dir.mkdir(exist_ok=True)

    ref_cds_fa = output / "reference_cps_cds.filtered.fa"
    ref_gene_order = extract_reference_cds(ref_gb, ref_cds_fa)

    mapping = blast_map(
        ref_cds_fa,
        pan_ref,
        threads=args.threads,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
    )
    if not mapping:
        raise RuntimeError("No CPS genes mapped from reference to Panaroo gene set")

    mapped_order = [(q, mapping[q]) for q in ref_gene_order if q in mapping]
    if not mapped_order:
        raise RuntimeError("BLAST mapped genes exist but none match reference order")

    # Preferred reconstruction source: gene_data.csv
    source = "gene_data"
    isolate_cols: List[str] = []
    gene_to_iso_seqlist: Dict[str, Dict[str, List[str]]] = {}

    if gene_data_csv.exists():
        isolates, g2i = read_gene_data(gene_data_csv)
        isolate_cols = sorted(isolates)
        gene_to_iso_seqlist = g2i
    else:
        source = "gene_presence_absence+combined_DNA_CDS"
        if not gpa_csv.exists() or not combined_cds.exists():
            raise FileNotFoundError(
                "Need either gene_data.csv OR both gene_presence_absence.csv and combined_DNA_CDS.fasta"
            )
        isolate_cols, gpa = read_gene_presence_absence(gpa_csv)
        seq_index = read_combined_dna_cds(combined_cds)

        g2i: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
        for gene, iso_map in gpa.items():
            for iso, ids in iso_map.items():
                seq = best_seq_for_ids(ids, seq_index)
                if seq:
                    g2i[gene][iso].append(seq)
        gene_to_iso_seqlist = g2i

    mapping_tsv = output / "cps_gene_mapping.tsv"
    with mapping_tsv.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["reference_gene", "panaroo_gene"])
        for rg, pg in mapped_order:
            w.writerow([rg, pg])

    usable_gene_alignments: List[Tuple[str, str, OrderedDict[str, str]]] = []
    all_isolates: Set[str] = set(isolate_cols)

    for ref_gene, pan_gene in mapped_order:
        iso_map = gene_to_iso_seqlist.get(pan_gene, {})
        per_isolate_seq: OrderedDict[str, str] = OrderedDict()
        for iso in isolate_cols:
            seq = best_seq_from_list(iso_map.get(iso, []))
            if seq:
                per_isolate_seq[iso] = seq

        if len(per_isolate_seq) < 2:
            continue

        raw_fa = raw_gene_dir / f"{pan_gene}.raw.fasta"
        write_fasta(raw_fa, {sanitize_isolate_name(k): v for k, v in per_isolate_seq.items()})

        aln_fa = cps_gene_align_dir / f"{pan_gene}.aln.fas"
        align_gene(raw_fa, aln_fa, args.threads)
        aln = read_alignment(aln_fa)
        usable_gene_alignments.append((ref_gene, pan_gene, aln))

    if not usable_gene_alignments:
        raise RuntimeError("No usable CPS genes could be reconstructed from selected Panaroo source")

    isolate_concat: Dict[str, List[str]] = {
        sanitize_isolate_name(iso): [] for iso in sorted(all_isolates)
    }

    for _, _, aln in usable_gene_alignments:
        aln_len = len(next(iter(aln.values())))
        for iso in isolate_concat:
            isolate_concat[iso].append(aln.get(iso, "N" * aln_len))

    concat_fasta = output / "cps_core_alignment.fasta"
    concat_records = OrderedDict((iso, "".join(parts)) for iso, parts in isolate_concat.items())
    write_fasta(concat_fasta, concat_records)

    for iso, seq in concat_records.items():
        write_fasta(cps_fastas_dir / f"{iso}_cps.fa", {f"{iso}_cps": seq})

    vcf_path = output / "cps_core_snps.vcf"
    run_cmd(["snp-sites", "-v", "-o", str(vcf_path), str(concat_fasta)])

    summary = output / "summary.tsv"
    with summary.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["reconstruction_source", source])
        w.writerow(["reference_genes_kept", len(ref_gene_order)])
        w.writerow(["mapped_genes", len(mapped_order)])
        w.writerow(["genes_with_alignments", len(usable_gene_alignments)])
        w.writerow(["isolates", len(concat_records)])

    print(f"Done. Results in: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
