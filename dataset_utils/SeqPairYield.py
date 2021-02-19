import numpy as np


class SeqPairYield:
    def __init__(self, thds, score):
        self.thd1 = thds[0]
        self.thd2 = thds[1]
        self.nupack_yield = score

    def as_list(self):
        single_cn = self.thd1.single_cn + self.thd2.single_cn
        single_mfe = self.thd1.single_mfe + self.thd2.single_mfe
        pair_cn = self.thd1.pair_cn + self.thd2.pair_cn
        pair_mfe = self.thd1.pair_mfe + self.thd2.pair_mfe
        label = "1" if self.nupack_yield > 0.2 else "0"
        gc_mean = np.mean([self.thd1.gc, self.thd2.gc])

        return [
            self.thd1.name,
            self.thd2.name,
            self.thd1.seq,
            self.thd2.seq,
            single_cn,
            single_mfe,
            pair_cn,
            pair_mfe,
            gc_mean,
            label,
            self.nupack_yield,
        ]
