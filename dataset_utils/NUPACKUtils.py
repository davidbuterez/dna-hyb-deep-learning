import subprocess
import shlex
from FileUtils import makedir, cd
from SeqPairYield import SeqPairYield


def complexes(path, temp=37.0):
    cmds = "complexes -T %.2f -material dna %s" % (temp, path)
    return (subprocess.run(shlex.split(cmds), stdout=subprocess.PIPE)).stdout


def concentrations(path):
    cmds = "concentrations %s" % (path,)
    return (subprocess.run(shlex.split(cmds), stdout=subprocess.PIPE)).stdout


def partition(path):
    cmds = "pfunc -material dna -multi %s" % (path,)
    return (subprocess.run(shlex.split(cmds), stdout=subprocess.PIPE)).stdout


def extract_mfe(mfe_file):
    mfe = []

    start_counting = False
    count = 0

    for line in mfe_file:
        if start_counting:
            count += 1

        if "composition" in line.strip():
            start_counting = True

        if count == 2:
            mfe.append(float(line.strip()))
            start_counting = False
            count = 0

    return mfe


def save_in_thermo(thd):
    seq = thd.seq
    name = "exp%s" % (seq,)
    makedir(name)

    with cd(name):
        with open("%s.in" % (name,), "w") as fh:
            print(1, file=fh)
            print(seq, file=fh)
            print(2, file=fh)
        with open("%s.con" % (name,), "w") as fh:
            print("1e-6", file=fh)
    return name


def nupack_thermo(complexname, thd):
    cn_score = 0.0
    mfe_score = 0.0

    cn_pair_score = 0.0
    mfe_pair_score = 0.0

    with cd(complexname):
        complexes(complexname)
        concentrations(complexname)
        with open("%s.eq" % (complexname,), "r") as eq_file:
            for line in eq_file:
                if not line.startswith("%"):
                    ln = line.strip().split("\t")
                    if ln[2] == "1":
                        cn_score = float(ln[-1]) * float("1e6")
                    elif ln[2] == "2":
                        cn_pair_score = float(ln[-1]) * float("1e6")

        with open("%s.ocx-mfe" % (complexname,), "r") as mfe_file:
            try:
                mfe = extract_mfe(mfe_file)
                mfe_score = mfe[0]
                mfe_pair_score = mfe[1]
            except:
                mfe_score = 0.0
                mfe_pair_score = 0.0

    thd.single_cn = cn_score
    thd.single_mfe = mfe_score
    thd.pair_cn = cn_pair_score
    thd.pair_mfe = mfe_pair_score


def nupack_thermo_analysis(thd):
    makedir("nupack_thermo2")
    with cd("nupack_thermo2"):
        complexname = save_in_thermo(thd)
        nupack_thermo(complexname, thd)
    return thd


def save_in_yield(thds):
    seq1 = thds[0].seq
    seq2 = thds[1].seq

    name = "exp%s%s" % (seq1, seq2)
    makedir(name)

    with cd(name):
        with open("%s.in" % (name,), "w") as fh:
            print(2, file=fh)
            print(seq1, file=fh)
            print(seq2, file=fh)
            print(2, file=fh)

        with open("%s.con" % (name,), "w") as fh:
            print("1e-6", file=fh)
            print("1e-6", file=fh)
    return name


def nupack_yield(complexname, t):
    with cd(complexname):
        complexes(complexname, temp=t)
        concentrations(complexname)
        with open("%s.eq" % (complexname,), "r") as eq_file:
            for line in eq_file:
                if not line.startswith("%"):
                    ln = line.split("\t")
                    if ln[2] == "1" and ln[3] == "1":
                        return float(ln[-2]) * float("1e6")
    return 0.0


def nupack_yield_analysis(thds, temp):
    makedir("nupack")
    with cd("nupack"):
        complexname = save_in_yield(thds)
        score = nupack_yield(complexname, temp)
        return SeqPairYield(thds, score)
