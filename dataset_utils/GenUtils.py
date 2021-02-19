import random

alphabet = set(["A", "C", "G", "T"])
ins_len_sample = [1, 2, 3, 4, 5]
del_len_sample = [1, 2, 3, 4]


def rev_comp(seq):
    rev = list(seq)[::-1]
    for i, ch in enumerate(rev):
        if ch == "A":
            rev[i] = "T"
        elif ch == "T":
            rev[i] = "A"
        elif ch == "C":
            rev[i] = "G"
        else:
            rev[i] = "C"
    return "".join(rev)


def seq(length):
    seq = []
    for i in range(length):
        if i > 1 and seq[i - 1] in ["G", "C"] and seq[i - 1] == seq[i - 2]:
            limited = alphabet.copy()
            limited = limited.difference(set([seq[i - 1]]))
            seq.append(random.sample(limited, 1)[0])
        elif (
            i > 2
            and seq[i - 1] in ["A", "T"]
            and seq[i - 1] == seq[i - 2]
            and seq[i - 2] == seq[i - 3]
        ):
            limited = alphabet.copy()
            limited = limited.difference(set([seq[i - 1]]))
            seq.append(random.sample(limited, 1)[0])
        else:
            seq.append(random.sample(alphabet, 1)[0])
    return "".join(seq)


def insertion(target, times):
    # MAX_TIMES = 3

    # Randomly choose how many insertions we do
    limit = random.randint(1, times)
    # print("Insert %d times." % (limit,))

    # Overlap keeps track of used indices
    overlap = set()

    # A copy of the sequence; to be modified
    alteration = list(target).copy()

    # Perform *limit* insertions (typically 1-3) per sequence
    for i in range(limit):
        indices = set(range(len(alteration) + 1))
        # Choose starting index and insertion length
        # Small insertions (typically 1-3) are more probable
        insert_index = random.sample(indices, 1)[0]
        insert_len = random.choice(ins_len_sample)

        # Check that a previous insertions has not occured at the same place
        while not overlap.isdisjoint(
            set(range(insert_index, insert_index + insert_len))
        ):
            insert_index = random.sample(indices, 1)[0]
            insert_len = random.choice(ins_len_sample)

        # Choose the nucleotides to be inserted
        insertion = [random.sample(alphabet, 1)[0] for _ in range(insert_len)]

        # Update the overlap set with the used indices
        overlap.update(
            list([i for i in range(insert_index, (insert_index + insert_len))])
        )

        # Debug info
        # print("Chosen insert index %d." % (insert_index,))
        # print("Overlap is: ", ", ".join(str(o) for o in overlap))
        # print("Chosen insertion %s of length %d." % ("".join(insertion),
        # insert_len))

        # Finally, insert the chosen nucleotides at the chosen position
        alteration[insert_index:insert_index] = insertion

    return "".join(alteration)


def deletion(target, times):
    # MAX_TIMES = 3

    # Randomly choose how many deletions we do
    limit = random.randint(1, times)

    # A copy of the sequence; to be modified
    alteration = list(target).copy()

    for i in range(limit):
        indices = set(range(len(alteration) + 1))
        # Choose starting index and deletion length
        # Small deletions (typically 1-3) are more probable
        delete_index = random.sample(indices, 1)[0]
        delete_len = random.choice(del_len_sample)

        if len(alteration) <= 1:
            return alteration

        # Check that we are not trying to delete out of bounds
        while delete_index + delete_len >= len(alteration):
            delete_index = random.sample(indices, 1)[0]
            delete_len = random.choice(del_len_sample)

        # Actual deletion
        del alteration[delete_index : (delete_index + delete_len)]

    return "".join(alteration)


def substitution(seq, times):
    # MAX_TIMES = 10

    # Randomly choose how many substitutions we do
    limit = random.randint(1, times)

    # A copy of the sequence; to be modified
    alteration = list(seq).copy()

    # Overlap keeps track of used indices
    overlap = set()

    for i in range(limit):
        indices = set(range(len(alteration)))
        # Choose substitution index
        subst_index = random.sample(indices, 1)[0]

        # Check we have not used this index before
        while subst_index in overlap:
            subst_index = random.sample(indices, 1)[0]

        # Add index to used indices set
        overlap.add(subst_index)

        alteration[subst_index] = random.sample(alphabet, 1)[0]

    return "".join(alteration)


def mutations(seq, fn, limits):
    muts = set()
    while len(muts) < limits[0]:
        muts.add(fn(seq, 1))

    while len(muts) < limits[1]:
        muts.add(fn(seq, 2))

    while len(muts) < limits[2]:
        muts.add(fn(seq, 3))
    return muts


def muts(seq, limits):
    subs = mutations(seq, substitution, limits)
    ins = mutations(seq, insertion, limits)
    dels = mutations(seq, deletion, limits)

    result = set()
    result.update(subs)
    result.update(ins)
    result.update(dels)
    return list(result)


def many_subs(seq, n):
    muts = set()

    times = [6, 7, 8, 9, 10, 11, 12]
    ins_times = [1, 2, 3]
    while len(muts) < n:
        muts.add(substitution(seq, random.choice(times)))
    while len(muts) < 2 * n:
        rnd = random.randint(0, 1)
        if rnd == 0:
            m = insertion(seq, random.choice(ins_times))
            if len(m) >= 18 and len(m) <= 26:
                muts.add(m)
        else:
            m = deletion(seq, random.choice(ins_times))
            if len(m) >= 18 and len(m) <= 26:
                muts.add(m)
    return muts
