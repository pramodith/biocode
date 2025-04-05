"""
Concept: RNA Hairpins and Secondary Structure
RNA strands can fold into complex shapes. One of the most common and biologically important motifs is the hairpin loop — a structure where a segment of RNA base-pairs with its reverse complement nearby, separated by a small unpaired loop in the middle.

A minimal hairpin loop consists of: - A stem of complementary base pairs (e.g., G–C, A–U) - A loop of 3 or more unpaired bases - The two stem regions are reverse complements of each other

Hairpin loops play roles in transcription regulation, stability, and even gene silencing (e.g., microRNAs).

Intuitive Example
Let's say we have the RNA strand:

GGGAAAUCCC

Left stem: GGG
Loop: AAAU
Right stem: CCC
If we reverse-complement the right stem: - CCC → GGG

And compare it to the left stem: - GGG vs GGG → ✅ match

This forms a valid hairpin.

Your Task
Write a function that takes an RNA strand and checks if it contains a valid hairpin loop. You'll be given: - The RNA string - The length of the stem (stem_len) - The length of the loop (loop_len)

A valid hairpin exists if: - The sequence contains stem_len bases - Followed by at least loop_len unpaired bases - Followed by stem_len bases that are the reverse complement of the initial stem
"""


def has_hairpin(seq: str, stem_len: int, loop_len: int) -> bool:
    # Your code here
    valid_stem = True
    m = {"G": "C", "A": "U", "C": "G", "U": "A"}

    for i in range(stem_len):
        if seq[i] != seq[-i - 1] and seq[i] == m[seq[-i - 1]]:
            continue
        else:
            valid_stem = False
            break
    l_c = 0
    for j in range(stem_len, len(seq)):
        if seq[j] == seq[-j - 1] or seq[j] != m[seq[-j - 1]]:
            l_c += 1

    return valid_stem and l_c == loop_len
