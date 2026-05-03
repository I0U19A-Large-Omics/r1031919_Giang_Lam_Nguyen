
#!/usr/bin/env python3



import argparse

import sqlite3

from pathlib import Path



import pandas as pd

import matplotlib.pyplot as plt





IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

COLORS = {

    "TLE66_N": "#0072B2",  # colorblind-safe blue

    "TLE66_T": "#D55E00",  # colorblind-safe orange

}





def query(db, sql):

    with sqlite3.connect(db) as con:

        return pd.read_sql_query(sql, con)





def plot_fig1(db, outfile):

    sql = """

    SELECT

        c.sample AS sample,

        e.annotation_impact AS impact,

        COUNT(DISTINCT s.chrom || ':' || s.pos) AS snp_count

    FROM SNP AS s

    JOIN Effect AS e ON s.snp_id = e.snp_id

    JOIN Call AS c ON s.snp_id = c.snp_id

    WHERE e.annotation_impact IN ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')

    GROUP BY c.sample, e.annotation_impact

    """



    df = query(db, sql)



    if df.empty:

        raise RuntimeError("Figure 1 query returned no rows.")



    # Ensure all impact/sample combinations exist, even if count is zero.

    samples = sorted(df["sample"].unique())

    full_index = pd.MultiIndex.from_product(

        [samples, IMPACT_ORDER],

        names=["sample", "impact"]

    )



    df = (

        df.set_index(["sample", "impact"])

          .reindex(full_index, fill_value=0)

          .reset_index()

    )



    fig, ax = plt.subplots(figsize=(8, 5))



    x = range(len(IMPACT_ORDER))

    width = 0.35



    for i, sample in enumerate(samples):

        sub = df[df["sample"] == sample]

        offsets = [j + (i - 0.5) * width for j in x]

        ax.bar(

            offsets,

            sub["snp_count"],

            width=width,

            label=sample,

            color=COLORS.get(sample, None),

        )



    ax.set_xticks(list(x))

    ax.set_xticklabels(IMPACT_ORDER)

    ax.set_xlabel("Predicted SNP impact category")

    ax.set_ylabel("Number of distinct SNP positions")

    ax.set_title("SNP impact severity per sample")

    ax.set_ylim(bottom=0)

    ax.legend(title="Sample")

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)



    fig.tight_layout()

    fig.savefig(outfile, format="svg")

    plt.close(fig)





def plot_fig2(db, outfile):

    # Horizontal bar chart: best for ranked categorical gene names, because labels stay readable.

    sql = """

    SELECT

        e.gene_name AS gene_name,

        COUNT(DISTINCT s.chrom || ':' || s.pos) AS snp_count

    FROM SNP AS s

    JOIN Effect AS e ON s.snp_id = e.snp_id

    WHERE e.annotation_impact IN ('HIGH', 'MODERATE')

      AND e.gene_name IS NOT NULL

      AND e.gene_name != ''

    GROUP BY e.gene_name

    ORDER BY snp_count DESC

    LIMIT 15

    """



    df = query(db, sql)



    if df.empty:

        raise RuntimeError("Figure 2 query returned no rows.")



    df = df.sort_values("snp_count", ascending=True)



    fig, ax = plt.subplots(figsize=(8, 6))



    ax.barh(df["gene_name"], df["snp_count"], color="#0072B2")

    ax.set_xlabel("Number of HIGH or MODERATE impact SNPs")

    ax.set_ylabel("Gene")

    ax.set_title("Genes most affected by predicted functional SNPs")

    ax.set_xlim(left=0)

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)



    fig.tight_layout()

    fig.savefig(outfile, format="svg")

    plt.close(fig)





def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--db", required=True)

    parser.add_argument("--fig1", required=True)

    parser.add_argument("--fig2", required=True)

    args = parser.parse_args()



    Path(args.fig1).parent.mkdir(parents=True, exist_ok=True)

    Path(args.fig2).parent.mkdir(parents=True, exist_ok=True)



    plot_fig1(args.db, args.fig1)

    plot_fig2(args.db, args.fig2)





if __name__ == "__main__":

    main()

