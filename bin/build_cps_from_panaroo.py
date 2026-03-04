#!/usr/bin/env python3
"""
Build CPS-focused outputs from Panaroo output and a serotype reference GenBank.

Workflow:
1) Extract reference CPS CDS sequences from GenBank, excluding pseudogenes/transposons.
2) BLAST these CDS against panaroo pan_genome_reference.fa to identify CPS genes robustly by sequence.
3) Subset Panaroo per-gene alignments for the matched genes.
4) Build per-isolate CPS FASTA (concatenated genes in reference order; missing genes become Ns).
5) Build concatenated CPS alignment + VCF (via snp-sites).
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature


EXCLUDE_TERMS = {
    "transposase",
    "transposon",
    "insertion sequence",
    "integrase",
    "recombinase",
}


def run_cmd(cmd: List[str], *, cwd: Path | None = None) -> None:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panaroo-dir", required=True, help="Path to Panaroo output directory")
    p.add_argument("--reference-gb", required=True, help="Reference serotype GenBank file")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--min-pident", type=float, default=90.0, help="Minimum BLAST percent identity")
    p.add_argument("--min-qcov", type=float, default=0.8, help="Minimum query coverage (0-1)")
    p.add_argument("--threads", type=int, default=4, help="Threads for BLAST")
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


def read_fasta_to_dict(path: Path) -> OrderedDict[str, str]:
    d: OrderedDict[str, str] = OrderedDict()
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = str(rec.seq)
    return d


def write_fasta(path: Path, seqs: Dict[str, str]) -> None:
    with path.open("w") as h:
        for name, seq in seqs.items():
            h.write(f">{name}\n{seq}\n")


def sanitize_isolate_name(name: str) -> str:
    # Keep behavior conservative and filesystem-friendly
    return name.replace("/", "_").replace(" ", "_")


def main() -> int:
    args = parse_args()

    panaroo_dir = Path(args.panaroo_dir)
    output = Path(args.output)
    ref_gb = Path(args.reference_gb)

    pan_ref = panaroo_dir / "pan_genome_reference.fa"
    aln_dir = panaroo_dir / "aligned_gene_sequences"

    if not pan_ref.exists():
        raise FileNotFoundError(f"Missing Panaroo reference FASTA: {pan_ref}")
    if not aln_dir.exists():
        raise FileNotFoundError(f"Missing Panaroo alignments directory: {aln_dir}")

    output.mkdir(parents=True, exist_ok=True)
    cps_gene_align_dir = output / "cps_gene_alignments"
    cps_gene_align_dir.mkdir(exist_ok=True)
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

    # preserve reference order, keeping only mapped genes
    mapped_order = [(q, mapping[q]) for q in ref_gene_order if q in mapping]

    presence_tsv = output / "cps_gene_mapping.tsv"
    with presence_tsv.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["reference_gene", "panaroo_gene"])
        for q, s in mapped_order:
            w.writerow([q, s])

    # Load alignments per mapped Panaroo gene
    gene_alignments: List[Tuple[str, str, OrderedDict[str, str]]] = []
    all_isolates = set()

    for ref_gene, pan_gene in mapped_order:
        aln_path = aln_dir / f"{pan_gene}.aln.fas"
        if not aln_path.exists():
            # keep pipeline going, but skip missing alignment files
            continue
        seqs = read_fasta_to_dict(aln_path)
        if not seqs:
            continue
        all_isolates.update(seqs.keys())
        gene_alignments.append((ref_gene, pan_gene, seqs))
        shutil.copy2(aln_path, cps_gene_align_dir / aln_path.name)

    if not gene_alignments:
        raise RuntimeError("No mapped CPS genes had usable Panaroo alignments")

    # Build per-isolate concatenated CPS
    isolate_concat: Dict[str, List[str]] = {iso: [] for iso in sorted(all_isolates)}

    for _, _, seqs in gene_alignments:
        aln_len = len(next(iter(seqs.values())))
        for iso in isolate_concat:
            isolate_concat[iso].append(seqs.get(iso, "N" * aln_len))

    concat_fasta = output / "cps_core_alignment.fasta"
    concat_records = OrderedDict(
        (sanitize_isolate_name(iso), "".join(parts)) for iso, parts in isolate_concat.items()
    )
    write_fasta(concat_fasta, concat_records)

    # Also write per-isolate cps FASTA
    for iso, seq in concat_records.items():
        write_fasta(cps_fastas_dir / f"{iso}_cps.fa", {f"{iso}_cps": seq})

    # VCF from concatenated alignment
    vcf_path = output / "cps_core_snps.vcf"
    run_cmd(["snp-sites", "-v", "-o", str(vcf_path), str(concat_fasta)])

    # Summary
    summary = output / "summary.tsv"
    with summary.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["reference_genes_kept", len(ref_gene_order)])
        w.writerow(["mapped_genes", len(mapped_order)])
        w.writerow(["genes_with_alignments", len(gene_alignments)])
        w.writerow(["isolates", len(concat_records)])

    print(f"Done. Results in: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
