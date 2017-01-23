# -*- coding: utf-8 -*-
"""
Code To Find Genes

@author: William Derksen

"""

import random
import math
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    if(nucleotide == 'A'):
        return 'T'
    if(nucleotide == 'T'):
        return 'A'
    if(nucleotide == 'C'):
        return 'G'
    if(nucleotide == 'G'):
        return 'C'
    pass


def get_reverse_complement(dna):
    newDna = ''
    for nucleo in dna:
        newNucleo = get_complement(nucleo)
        newDna = newDna + newNucleo
    return newDna


def rest_of_ORF(dna):
    active = False
    i = 0
    ORFLeft = ''
    """while i < len(dna)-2:
        nucleo = dna[i] + dna[i+1] + dna[i+2]
        if nucleo == 'TAG' or nucleo == 'TAA' or nucleo == 'TGA':
            return ORFLeft
        ORFLeft = ORFLeft + dna[i]
        i+=1
    return ORFLeft"""
    while i < len(dna):
        nucleo = dna[i:i+3]
        if active == True or nucleo == 'ATG':
            active = True
            if nucleo == 'TAG' or nucleo == 'TAA' or nucleo == 'TGA':
                active = False
                return [ORFLeft,i]
            ORFLeft = ORFLeft + nucleo
        i+=3
    if ORFLeft != '':
        return [ORFLeft,i]
    else:
        return [None, i]


print(rest_of_ORF('GCATGTATGCAG'))


def find_all_ORFs_oneframe(dna):
    ORFs = []
    i=0
    while i<len(dna):
        newOrf = rest_of_ORF(dna[i:])
        if newOrf[0] == None:
            return ORFs
        i+=newOrf[1]+3
        ORFs.append(newOrf[0])
    return ORFs
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    pass


print(find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC"))

def find_all_ORFs(dna):
    ORFs = []
    i = 0
    while i < 3:
        ORFs = ORFs + find_all_ORFs_oneframe(dna[i:])
        i+=1
    return ORFs
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    pass


print(find_all_ORFs("ATGCATGAATGTAG"))

def find_all_ORFs_both_strands(dna):
    otherRevDna = get_reverse_complement(dna[::-1])
    dnaORF = find_all_ORFs(dna)
    otherRevDnaOrf = find_all_ORFs(otherRevDna)
    return dnaORF + otherRevDnaOrf
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass

print(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))


def longest_ORF(dna):
    Orfs = find_all_ORFs_both_strands(dna)
    lengs = []
    for Orf in Orfs:
        lengs.append(len(Orf))
    maxlen = max(lengs)
    ind = lengs.index(maxlen) #IMPORTANT NOTE, ONLY RETURNS THE FIRST MAX FOUND, NOT ALL MAXS"""
    return Orfs[ind]

    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass

print(longest_ORF("ATGCGAATGTAGCATCAAA"))

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


print(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)
