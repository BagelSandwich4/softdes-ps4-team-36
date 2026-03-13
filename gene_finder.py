"""
Library for finding potential genes in a strand of DNA.
"""

import random
from helpers import amino_acid
from helpers import load_fasta_file


def get_complement(nucleotide):
    """
    Finds the complementary nucleotide to an input (ex. "A"->"T").

    ARGS:
        nucleotide, string, single capitol letter representing the nucleotide
        you want to find the completmentary version of (ex. "A")

    Returns:
        a string with a single capitol letter representing the complement
        of the input nucleotide (ex. "T").
    """
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "C":
        return "G"
    if nucleotide == "G":
        return "C"
    return ""


def get_reverse_complement(strand):
    """
    Finds the reverse complement strand from any given DNA strand.

    ARGS:
        strand, string, string consisting only of the following
        characters: "A", "T", "G", "C" representing a strand of DNA.

    Returns:
        string, reverse complement of input strand.
    """
    strand_list = list(strand)
    reverse_strand = []
    for char in strand_list:
        reverse_strand.append(get_complement(char))
    reverse_strand.reverse()
    return "".join(reverse_strand)


def rest_of_orf(strand):
    """
    Takes a strand of DNA that begins with a start codon and returns the
    sequence of nucleotides representing the rest of the ORF up to the
    next stop codon (non-inclusive). If an ORF is not provided, this will
    return nothing.

    Args:
        strand, string representing DNA as a sequence of nucleotides.

    Returns:
        string of the ORF from the start to the next stop codon.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(0, len(strand) - 2, 3):
        codon = strand[i : i + 3]
        if codon in stop_codons:
            return strand[:i]
    return ""


def find_all_orfs_one_frame(strand):
    """
    Takes string strand and returns a list representing all
    in-frame ORFs.

    Args:
        strand, string representing DNA as a sequence of nucleotides.

    Returns:
        list, each element is a string representing an in-frame
        ORF of the givne strand.
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
    Finds all ORFs of a given strand, returns all ORFs found in
    that strand including not just in-frame ORFs.

    Args:
        strand, string representing a strand of DNA

    Returns:
        orfs, list of strings representing all ORFs found in the
        strand of DNA, including those out-of-frame.
    """
    orfs = []

    orfs.extend(find_all_orfs_one_frame(strand))
    orfs.extend(find_all_orfs_one_frame(strand[1:]))
    orfs.extend(find_all_orfs_one_frame(strand[2:]))

    return orfs


def find_all_orfs_both_strands(strand):
    """
    takes a strand of DNA and returns a list of strings representing
    all ORFs found in the strand or its reverse complement.

    Args:
        strand, string representing a sequence of DNA.

    Returns:
        orfs: a list of strings representing all ORFs found in the
        given DNA strand as well is its reverse complement strand.
    """
    orfs = []
    reverse = get_reverse_complement(strand)
    orfs.extend(find_all_orfs(reverse))
    orfs.extend(find_all_orfs(strand))
    return orfs


def find_longest_orf(strand):
    """
    Returns the longest ORF found in a given strand as well as in its
    reverse complement strand.

    Args:
        strand, string representing a sequence of DNA.

    Returns:
        string representing the longest ORF found in the DNA sequence.
    """
    orfs = find_all_orfs_both_strands(strand)
    if not orfs:
        return ""
    return max(orfs, key=len)


def noncoding_orf_threshold(strand, num_trials):
    """
    Returns the lowest length ORF in randomly shuffled DNA strands
    over multiple random strands.

    Args:
        strand, string representing a strand of DNA

        num_trials, int representing how many random strains are
        being tested

    Returns:
        int representing the ratio of the minimum ORF length to the longest
        ORF length found by randomly shuffling the input DNA strand.
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
    Converts a given ORF into a sequence of amino acids.

    Args:
        orf, string representing a strand of DNA that is an ORF.

    Returns:
        result, string representing an amino acid sequence.
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
    Finds potential protein-coding genes in a string of nucleotides
    from a file in a NIH database, then returns a list of all amino
    acid sequences in the longest ORF strands from 1500 randomly
    shuffled trials of the target strand, found in the database at
    path

    Args:
        path, string representing the location of a file in the NIH
        databse

    Returns:
        genes, list of all amino acid sequences
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
