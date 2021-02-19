import pandas as pd
import parasail
from SeqThermoData import SeqThermoData
from NUPACKUtils import nupack_thermo_analysis
from pathos.multiprocessing import ProcessPool
from Bio.Seq import Seq

df = pd.read_csv("ds.csv")

seq_dict = {}
for index, row in df.iterrows():
    s1 = row.Seq1
    s2 = row.Seq2
    if s1 not in seq_dict:
        seq_dict[s1] = SeqThermoData(None, s1)
    if s2 not in seq_dict:
        seq_dict[s2] = SeqThermoData(None, s2)


with ProcessPool() as p:
    thermo_data = p.map(nupack_thermo_analysis, seq_dict.values())

for thd in thermo_data:
    seq_dict[thd.seq] = thd


dataset = []
for i, row in df.iterrows():
    thd1 = seq_dict[row.Seq1]
    thd2 = seq_dict[row.Seq2]

    label = 1 if row.Yield > 0.2 else 0

    s1 = Seq(thd1.seq)
    s1 = s1.reverse_complement()
    aln = parasail.sg_scan_16(str(s1), row.Seq2, 5, 2, parasail.dnafull)
    dataset.append(
        [
            thd1.seq,
            thd2.seq,
            aln.score,
            thd1.gc,
            thd2.gc,
            thd1.single_cn,
            thd2.single_cn,
            thd1.pair_cn,
            thd2.pair_cn,
            thd1.single_mfe,
            thd2.single_mfe,
            row.Yield,
            label,
        ]
    )


dff = pd.DataFrame(
    dataset,
    columns=[
        "Seq1",
        "Seq2",
        "Aln",
        "GC1",
        "GC2",
        "SingleCon1",
        "SingleCon2",
        "PairCon1",
        "PairCon2",
        "SingleMFE1",
        "SingleMFE2",
        "Yield",
        "Label",
    ],
)


dff.to_csv("AlignSet3.csv", index=None)

