from Bio.Seq import Seq
from Bio.SeqUtils import GC


class SeqThermoData:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.gc = GC(self.seq)
        self.single_cn = 0.0
        self.single_mfe = 0.0
        self.pair_cn = 0.0
        self.pair_mfe = 0.0

    def __str__(self):
        return "%s\t%.8f\t%.8f\t%.8f\t%.8f" % (
            self.seq,
            self.single_cn,
            self.single_mfe,
            self.pair_cn,
            self.pair_mfe,
        )

    def get_rev_comp(self):
        s = Seq(self.seq)
        return s.reverse_complement()
