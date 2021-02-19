from GenUtils import seq, many_subs, rev_comp
from NUPACKUtils import nupack_yield_analysis
from SeqThermoData import SeqThermoData
from pathos.multiprocessing import ProcessPool
import random


def gen(i):
    lens = [18, 19, 20, 21, 22, 23, 24, 25, 26]

    sq = seq(random.choice(lens))
    rsq = rev_comp(sq)
    mts = many_subs(rsq, 100)

    pairs = []
    for m in mts:
        if len(m) >= 18 and len(m) <= 26:
            pairs.append((SeqThermoData(None, sq), SeqThermoData(None, m)))

    for pr in pairs:
        nupack_yield_analysis(pr, 57.0)


with ProcessPool() as p:
    result = p.map(gen, range(1000))
