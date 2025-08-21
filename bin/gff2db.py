#!/usr/bin/env python3
import os
import json
import logging
import argparse
from typing import Dict, List, Tuple

import screed

# ------------------------------- #
# Logging
# ------------------------------- #
def setup_logging(level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("fungalAMR.gff2db")
    if not logger.handlers:
        logger.setLevel(getattr(logging, level.upper(), logging.INFO))
        h = logging.StreamHandler()
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                                datefmt="%Y-%m-%d %H:%M:%S")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger

LOGGER = setup_logging()


# ------------------------------- #
# CLI
# ------------------------------- #
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create fungalAMR database from a GFF file and associated assembly FASTA."
    )
    p.add_argument("--gff", required=True, help="Path to GFF3 file.")
    p.add_argument("--fasta", required=True, help="Path to assembly FASTA file.")
    p.add_argument("--outdir", default=".", help="Output directory.")
    p.add_argument("--gene", help="Gene target (e.g., FKS1). If omitted, output all genes.")
    p.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR).")
    return p.parse_args()


# ------------------------------- #
# Core helpers
# ------------------------------- #
def load_genome_sequences(fasta_path: str) -> Dict[str, str]:
    """Load genome sequences keyed by the first token of each FASTA header."""
    LOGGER.info(f"Loading genome sequences from: {fasta_path}")
    genome: Dict[str, str] = {}
    n = 0
    with screed.open(fasta_path) as fh:
        for rec in fh:
            key = rec["name"].split()[0]
            genome[key] = rec["sequence"]
            n += 1
    LOGGER.info(f"Loaded {n} sequences from FASTA.")
    return genome


def parse_gff(gff_path: str) -> Dict[str, Dict[str, dict]]:
    """Parse GFF3 into feature buckets we care about: gene, mRNA, exon, CDS."""
    LOGGER.info(f"Parsing GFF: {gff_path}")
    data: Dict[str, Dict[str, dict]] = {"gene": {}, "mRNA": {}, "exon": {}, "CDS": {}}
    kept, skipped, total = 0, 0, 0

    with open(gff_path, "r", encoding="utf-8") as f:
        for line in f:
            total += 1
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                skipped += 1
                continue

            acc, source, ftype, start, end, score, strand, phase, attr = parts
            if ftype not in data:
                continue

            # Parse attributes
            attrs_lc: Dict[str, str] = {}
            for p in attr.split(";"):
                if not p:
                    continue
                k_v = p.split("=", 1)
                if len(k_v) != 2:
                    continue
                k, v = k_v
                attrs_lc[k.lower()] = v

            rec_id = attrs_lc.get("id")
            if not rec_id:
                skipped += 1
                continue

            rec = {
                "acc": acc,
                "ftype": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phase": None if phase == "." else phase,
                # minimal wiring attributes (lowercase keys)
                "id": rec_id,
                "parent": attrs_lc.get("parent"),
                "gene": attrs_lc.get("gene"),
                "product": attrs_lc.get("product"),
                "name": attrs_lc.get("name"),
            }
            data[ftype][rec_id] = rec
            kept += 1

    LOGGER.info(f"GFF parsed. Kept: {kept}, Skipped: {skipped}, Total lines: {total}")
    LOGGER.debug(
        "Counts by type: gene=%d, mRNA=%d, exon=%d, CDS=%d",
        len(data["gene"]), len(data["mRNA"]), len(data["exon"]), len(data["CDS"])
    )
    return data


def build_rna_map(data: Dict[str, Dict[str, dict]]) -> Dict[str, List[dict]]:
    """Build mRNA -> [mRNA, gene, exons, CDS] map."""
    LOGGER.info("Building mRNA-to-features map.")
    rna_map: Dict[str, List[dict]] = {}

    # CDS/exon features onto each parent mRNA
    for feat in data["CDS"].values():
        for parent in (feat.get("parent") or "").split(","):
            parent = parent.strip()
            if parent:
                rna_map.setdefault(parent, []).append(feat)

    for feat in data["exon"].values():
        for parent in (feat.get("parent") or "").split(","):
            parent = parent.strip()
            if parent:
                rna_map.setdefault(parent, []).append(feat)

    # Include the mRNA itself and its parent gene
    for mrna_id, mrna in data["mRNA"].items():
        rna_map.setdefault(mrna_id, []).append(mrna)
        gene_id = mrna.get("parent")
        if gene_id and gene_id in data["gene"]:
            rna_map[mrna_id].append({**data["gene"][gene_id]})

    LOGGER.info(f"mRNA map built for {len(rna_map)} transcripts.")
    return rna_map


def _normalize_coords(
    f_start: int, f_end: int, rec_start: int, rec_end: int, strand: str
) -> Tuple[int, int]:
    """Return (rel_start, rel_end) normalized to the transcript span."""
    if strand == "-":
        rel_start = rec_end - f_end + 1
        rel_end = rec_end - f_start + 1
    else:
        rel_start = f_start - rec_start + 1
        rel_end = f_end - rec_start + 1
    return rel_start, rel_end


def make_gene_objects(
    rna_map: Dict[str, List[dict]],
    data: Dict[str, Dict[str, dict]],
    genome: Dict[str, str],
    target_gene: str | None,
) -> List[dict]:
    """Construct gene objects per transcript; filter by target if provided."""
    LOGGER.info("Constructing gene objects.")
    genes: List[dict] = []
    built, filtered_missing_acc, filtered_no_mrna = 0, 0, 0

    for mrna_id, features in rna_map.items():
        mrna = data["mRNA"].get(mrna_id)
        if not mrna:
            filtered_no_mrna += 1
            continue

        gene_id = mrna.get("parent")
        gene_rec = data["gene"].get(gene_id, {})

        # Span for this transcript bundle
        rec_start = min(int(f["start"]) for f in features)
        rec_end = max(int(f["end"]) for f in features)
        strand = gene_rec.get("strand") or mrna.get("strand")
        acc = gene_rec.get("acc") or mrna.get("acc")

        seq = genome.get(acc)
        if not seq:
            filtered_missing_acc += 1
            continue

        # Normalize features and strip attributes
        norm_features: List[dict] = []
        for f in features:
            rs, re = _normalize_coords(int(f["start"]), int(f["end"]), rec_start, rec_end, strand)
            out_feat = {"type": f["ftype"], "coords": [rs, re]}
            if f.get("phase") not in (None, ".", ""):
                out_feat["phase"] = f["phase"]
            norm_features.append(out_feat)

        # Extract genomic sequence for the span (end inclusive in GFF)
        subseq = seq[rec_start - 1 : rec_end].upper()

        gene_out = {
            "acc": acc,
            "gene": gene_rec.get("gene"),
            "product": (gene_rec.get("product") or "None"),
            "strand": strand,
            "length": len(subseq),
            "sequence": subseq,
            "features": norm_features,
            "targets": [{"name": "target_1", "coords": []}],
        }

        # Filter by target (if provided)
        if target_gene:
            if (gene_out["gene"] or "").upper() == target_gene.upper():
                genes.append(gene_out)
                built += 1
        else:
            genes.append(gene_out)
            built += 1

    LOGGER.info(f"Constructed {built} gene objects.")
    if filtered_missing_acc:
        LOGGER.warning(f"Skipped {filtered_missing_acc} transcripts due to missing FASTA accession.")
    if filtered_no_mrna:
        LOGGER.warning(f"Skipped {filtered_no_mrna} entries lacking mRNA records.")
    return genes


def write_output(genes: List[dict], outdir: str, filename: str = "genes.json") -> str:
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, filename)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(genes, f, indent=2)
    LOGGER.info(f"Wrote {len(genes)} records to {out_path}")
    return out_path


# ------------------------------- #
# Main
# ------------------------------- #
def main() -> None:
    args = parse_args()
    # Update logging level now that we have args
    LOGGER.setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    LOGGER.debug(f"Arguments: {vars(args)}")

    genome = load_genome_sequences(args.fasta)
    data = parse_gff(args.gff)
    rna_map = build_rna_map(data)
    genes = make_gene_objects(rna_map, data, genome, args.gene)
    write_output(genes, args.outdir)


if __name__ == "__main__":
    main()
