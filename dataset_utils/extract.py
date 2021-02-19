import os
import pandas as pd
import numpy as np
from FileUtils import cd
from pathos.multiprocessing import ProcessPool


def process_dirs(dir_lst):
    score = 0.0
    records = []
    for dirname in dir_lst:
        with cd(dirname):
            seqs = []
            with open("%s.in" % (dirname,), "r") as infile:
                for line in infile:
                    read = line.strip().split()
                    if read[0] != "2":
                        seqs.append(read[0])
            try:
                with open("%s.eq" % (dirname,), "r") as eqfile:
                    for line in eqfile:
                        if not line.startswith("%"):
                            ln = line.split("\t")
                            if ln[2] == "1" and ln[3] == "1":
                                score = float(ln[-2]) * float("1e6")
            except:
                continue

            label = 1 if score > 0.2 else 0
            if score <= 1.0:
                try:
                    records.append([seqs[0], seqs[1], label, score])
                except:
                    print("DIRNAME %s" % (dirname,))
                    continue
    return records


dirs = os.listdir("nupack")
splits = np.array_split(dirs, len(dirs) / 200)

with cd("nupack"):
    with ProcessPool() as p:
        result = p.map(process_dirs, splits)


flat_list = [item for sublist in result for item in sublist]

headers = ["Seq1", "Seq2", "Label", "Yield"]
df = pd.DataFrame(flat_list, columns=headers)
df.to_csv("records_o_s.csv", index=None)
