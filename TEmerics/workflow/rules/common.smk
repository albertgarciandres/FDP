import pandas as pd

from pathlib import Path

###############################################################################
### Functions
###############################################################################

def get_sample(column_id: str, sample_id: int = None) -> str:
    """Get relevant per sample information."""
    if sample_id:
        return str(samples_table[column_id][samples_table.index == sample_id].iloc[0])
    else:
        return str(samples_table[column_id].iloc[0])


def get_trinity_input(df, sample_dir: Path, out_dir: Path) -> list(str):
    """Generate trinity input file from samples table."""

    strs = set(df["strain"])
    tiss = set(df["tissue"])

    smp = dict()
    conditions = list()
    
    for t in tiss:
        for s in strs:
            smp[f"{tissue}_{strain}"] = df.loc[(df["tissue"] == t) & (df["strain"] == s), "sample"].to_list()

    for cond in smp.keys():
        entries = smp[cond]

        if len(entries) > 0:
            conditions.append(cond)

            with open(f"{out_dir}/trinity_{condition}.tsv", "w") as out_file:

                for entry in smp[cond]:
                    path1 = f"{sample_dir}/{entry}_reads_trimmed_1.fq.gz"
                    path2 = f"{sample_dir}/{entry}_reads_trimmed_2.fq.gz"
                    out_file.write('\t'.join(cond, entry, path1, path2))

    return conditions

