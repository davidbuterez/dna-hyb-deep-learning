from GenUtils import seq, muts, rev_comp
from NUPACKUtils import nupack_yield_analysis
from SeqThermoData import SeqThermoData
from pathos.multiprocessing import ProcessPool
import random


def gen(i):
    lens = [18, 19, 20, 21, 22, 23, 24, 25, 26]

    sq = seq(random.choice(lens))
    mts = muts(sq, [10, 20, 30])

    pairs = []
    for m in mts:
        if len(m) >= 18 and len(m) <= 26:
            pairs.append((SeqThermoData(None, sq), SeqThermoData(None, m)))

    for pr in pairs:
        nupack_yield_analysis(pr, 57.0)


with ProcessPool() as p:
    result = p.map(gen, range(5000))
