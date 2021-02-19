from SeqThermoData import SeqThermoData
from Bio.Seq import Seq


def write_fasta(seqs, name):
    id_to_seq = {"%s%d" % (name, i): sq for i, sq in enumerate(list(seqs))}

    filename = "%s.fasta" % (name,)
    # Write seqs to FASTA file
    with open(filename, "w") as seqsfile:
        for name, sq in id_to_seq.items():
            print(">%s\n%s" % (name, sq), file=seqsfile)

    return filename


def read_fasta(filename):
    thds = []
    with open(filename, "r") as fh:
        while True:
            id_line = fh.readline().strip()
            seq_line = fh.readline().strip()
            if not id_line:
                break
            thds.append(SeqThermoData(id_line[1:], seq_line))

    return thds


def rev_comp_fasta(filename):
    seqs = []
    with open(filename, "r") as fh:
        while True:
            id_line = fh.readline().strip()
            seq_line = fh.readline().strip()
            if not id_line:
                break
            seqs.append(seq_line)

    rev_seqs = [Seq(sq).reverse_complement() for sq in seqs]
    write_fasta(rev_seqs, "rev")
