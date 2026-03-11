"""
Library for finding potential genes in a strand of DNA.
"""

from helpers import amino_acid
from helpers import load_fasta_file
import random


def get_complement(nucleotide):
    """
    Your docstring goes here.
    """
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "C":
        return "G"
    if nucleotide == "G":
        return "C"


def get_reverse_complement(strand):
    """
    Your docstring goes here.
    """
    strand_list = list(strand)
    reverse_strand = []
    for char in strand_list:
        reverse_strand.append(get_complement(char))
    reverse_strand.reverse()
    return "".join(reverse_strand)


def rest_of_orf(strand):
    """
    Your docstring goes here.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(0, len(strand) - 2, 3):
        codon = strand[i : i + 3]
        if codon in stop_codons:
            return strand[:i]
    return ""


def find_all_orfs_one_frame(strand):
    """
    Your docstring goes here.
    """
    orfs = []
    i = 0
    while i < len(strand) - 2:
        codon = strand[i : i + 3]
        if codon == "ATG":
            orf = rest_of_orf(strand[i:])
            if orf != "":
                orfs.append(orf)
                i += len(orf)
                continue
        i += 3
    return orfs


def find_all_orfs(strand):
    """
    Your docstring goes here.
    """
    orfs = []

    orfs.extend(find_all_orfs_one_frame(strand))
    orfs.extend(find_all_orfs_one_frame(strand[1:]))
    orfs.extend(find_all_orfs_one_frame(strand[2:]))

    return orfs


def find_all_orfs_both_strands(strand):
    """
    Your docstring goes here.
    """
    orfs = []
    reverse = get_reverse_complement(strand)
    orfs.extend(find_all_orfs(reverse))
    orfs.extend(find_all_orfs(strand))
    return orfs


def find_longest_orf(strand):
    """
    Your docstring goes here.
    """
    orfs = find_all_orfs_both_strands(strand)
    if not orfs:
        return ""
    return max(orfs, key=len)


def noncoding_orf_threshold(strand, num_trials):
    """
    Your docstring goes here.
    """
    min_length = float("inf")
    for _ in range(num_trials):
        strand_list = list(strand)
        random.shuffle(strand_list)
        shuffled = "".join(strand_list)
        longest = find_longest_orf(shuffled)
        if longest != "":
            min_length = min(min_length, len(longest))
    return min_length


def encode_amino_acids(orf):
    """
    Your docstring goes here.
    """
    result = ""
    for i in range(0, len(orf), 3):
        codon = orf[i : i + 3]
        if len(codon) < 3:
            break
        aa = amino_acid(codon)
        if aa == "*":
            break
        result += aa
    return result


def find_genes(path):
    """
    Your docstring goes here.
    """
    strand = load_fasta_file(path)
    threshold = noncoding_orf_threshold(strand, 1500)
    orfs = find_all_orfs_both_strands(strand)
    genes = []
    for orf in orfs:
        if len(orf) > threshold:
            protein = encode_amino_acids(orf)
            genes.append(protein)
    return genes


# DON'T ADD ANYTHING ELSE TO THIS FILE. IF YOU DO, YOUR CODE MAY BE MARKED AS
# FAILING UNIT TESTS. IF YOU NEED TO TEST YOUR CODE, DO IT IN A JUPYTER
# NOTEBOOK.
