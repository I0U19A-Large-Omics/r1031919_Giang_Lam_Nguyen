#!/usr/bin/env python3
"""
Score all SNP coordinates from an annotated VCF with fake_enformer.

Parallelism approach:
This script uses asyncio with fake_enformer.async_predict. The fake Enformer
calls are slow per coordinate, so we launch many predictions concurrently while
limiting concurrency with a semaphore. This avoids a slow sequential loop and
keeps the output reproducible by sorting results back into VCF order.
"""

import argparse
import asyncio
import inspect
import time
from pathlib import Path

import fake_enformer


def parse_snps_from_vcf(vcf_path: str) -> list[str]:
    coords = []

    with open(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom, pos, _id, ref, alt = fields[:5]

            for alt_allele in alt.split(","):
                if len(ref) == 1 and len(alt_allele) == 1 and alt_allele not in {".", "*"}:
                    coords.append(f"hg38:{chrom}:{pos}:{ref}:{alt_allele}")
                    break

    return coords


async def score_one(coord: str, semaphore: asyncio.Semaphore) -> tuple[str, float]:
    async with semaphore:
        result = fake_enformer.async_predict(coord)

        if inspect.isawaitable(result):
            result = await result

        return coord, float(result)


async def score_all(coords: list[str], concurrency: int) -> list[tuple[str, float]]:
    semaphore = asyncio.Semaphore(concurrency)
    tasks = [asyncio.create_task(score_one(coord, semaphore)) for coord in coords]

    results = []
    for i, task in enumerate(asyncio.as_completed(tasks), start=1):
        results.append(await task)
        if i % 10 == 0 or i == len(tasks):
            print(f"Completed {i}/{len(tasks)} SNPs", flush=True)

    order = {coord: i for i, coord in enumerate(coords)}
    results.sort(key=lambda item: order[item[0]])
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Input annotated VCF")
    parser.add_argument("-o", "--output", default="enformer_scores.tsv")
    parser.add_argument("-t", "--threads", type=int, default=16)
    args = parser.parse_args()

    start = time.perf_counter()

    coords = parse_snps_from_vcf(args.vcf)
    if not coords:
        raise SystemExit(f"No SNPs found in {args.vcf}")

    results = asyncio.run(score_all(coords, args.threads))

    out = Path(args.output)
    with out.open("wt", encoding="utf-8") as handle:
        handle.write("coordinate\tscore\n")
        for coord, score in results:
            handle.write(f"{coord}\t{score:.6g}\n")

    elapsed = time.perf_counter() - start
    print(f"VCF: {args.vcf}")
    print(f"SNPs scored: {len(results)}")
    print(f"Concurrency: {args.threads}")
    print(f"Output: {out}")
    print(f"Elapsed_seconds: {elapsed:.2f}")


if __name__ == "__main__":
    main()
