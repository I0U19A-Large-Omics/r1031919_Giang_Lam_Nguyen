#!/usr/bin/env python3
import argparse
import os
import sqlite3
import pandas as pd
import vcfpy

ANN_FIELDS = [
    "allele","annotation","annotation_impact","gene_name","gene_id",
    "feature_type","feature_id","transcript_biotype","rank_total",
    "hgvs_c","hgvs_p","cdna_pos_len","cds_pos_len","aa_pos_len",
    "distance","errors_warnings_info",
]

def parse_ann_entry(entry):
    parts = entry.split("|")
    if len(parts) < 16:
        parts.extend([""] * (16 - len(parts)))
    elif len(parts) > 16:
        parts = parts[:15] + ["|".join(parts[15:])]
    return dict(zip(ANN_FIELDS, parts))

def get_gt(call):
    if call.data is None:
        return None
    gt = call.data.get("GT")
    return None if gt is None else str(gt)

def flush_chunk(conn, snp_rows, effect_rows, call_rows, first_chunk):
    if not snp_rows:
        return first_chunk
    mode = "replace" if first_chunk else "append"
    pd.DataFrame(snp_rows).to_sql("SNP", conn, if_exists=mode, index=False)
    pd.DataFrame(effect_rows).to_sql("Effect", conn, if_exists=mode, index=False)
    pd.DataFrame(call_rows).to_sql("Call", conn, if_exists=mode, index=False)
    return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--chunk-size", type=int, default=500)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.db), exist_ok=True)

    reader = vcfpy.Reader.from_path(args.vcf)

    snp_rows, effect_rows, call_rows = [], [], []
    snp_id = 1
    effect_id = 1
    call_id = 1
    first_chunk = True

    with sqlite3.connect(args.db) as conn:
        for record in reader:
            alt = ",".join(str(a.value) for a in record.ALT)
            filt = None
            if record.FILTER:
                filt = ";".join(record.FILTER) if isinstance(record.FILTER, list) else str(record.FILTER)

            snp_rows.append({
                "snp_id": snp_id,
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": alt,
                "qual": record.QUAL,
                "filter": filt,
            })

            ann_entries = record.INFO.get("ANN", [])
            if isinstance(ann_entries, str):
                ann_entries = [ann_entries]

            for ann in ann_entries:
                parsed = parse_ann_entry(ann)
                effect_rows.append({
                    "effect_id": effect_id,
                    "snp_id": snp_id,
                    "allele": parsed["allele"],
                    "annotation": parsed["annotation"],
                    "annotation_impact": parsed["annotation_impact"],
                    "gene_name": parsed["gene_name"],
                    "gene_id": parsed["gene_id"],
                    "feature_type": parsed["feature_type"],
                    "feature_id": parsed["feature_id"],
                    "transcript_biotype": parsed["transcript_biotype"],
                    "rank_total": parsed["rank_total"],
                    "hgvs_c": parsed["hgvs_c"],
                    "hgvs_p": parsed["hgvs_p"],
                    "cdna_pos_len": parsed["cdna_pos_len"],
                    "cds_pos_len": parsed["cds_pos_len"],
                    "aa_pos_len": parsed["aa_pos_len"],
                    "distance": parsed["distance"],
                    "errors_warnings_info": parsed["errors_warnings_info"],
                })
                effect_id += 1

            for call in record.calls:
                call_rows.append({
                    "call_id": call_id,
                    "snp_id": snp_id,
                    "sample": call.sample,
                    "gt": get_gt(call),
                })
                call_id += 1

            snp_id += 1

            if len(snp_rows) >= args.chunk_size:
                first_chunk = flush_chunk(conn, snp_rows, effect_rows, call_rows, first_chunk)
                snp_rows, effect_rows, call_rows = [], [], []

        first_chunk = flush_chunk(conn, snp_rows, effect_rows, call_rows, first_chunk)

        conn.execute("CREATE INDEX IF NOT EXISTS idx_snp_snp_id ON SNP(snp_id);")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_effect_snp_id ON Effect(snp_id);")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_call_snp_id ON Call(snp_id);")
        conn.commit()

if __name__ == "__main__":
    main()
