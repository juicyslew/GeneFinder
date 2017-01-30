# -*- coding: utf-8 -*-
"""
Code To Find Genes

@author: William Derksen

"""

import random
import math
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq, load_nitrogenase_seq, load_metagenome
nitrogenase = load_nitrogenase_seq()
metagenome = load_metagenome()

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if(nucleotide == 'A'):
        return 'T'
    if(nucleotide == 'T'):
        return 'A'
    if(nucleotide == 'C'):
        return 'G'
    if(nucleotide == 'G'):
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    newDna = ''
    for nucleo in dna:
        newNucleo = get_complement(nucleo)
        newDna = newDna + newNucleo
    return newDna[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
        ##This function works differently from the intended ORF
    >>> rest_of_ORF("ATGTGAA")
    ['ATG', 3]
    >>> rest_of_ORF("ATGAGATAGG")
    ['ATGAGA', 6]
    """
    active = False
    i = 0
    ORFLeft = '' #Instantiate ORF
    while i < len(dna):
        nucleo = dna[i:i+3] #check groups of 3 letters
        if active == True or nucleo == 'ATG':  #if not active then only activate when "ATG" is found, after that stay active
            active = True
            if nucleo == 'TAG' or nucleo == 'TAA' or nucleo == 'TGA': #stop ORF
                active = False
                return [ORFLeft,i] # Return ORF
            ORFLeft = ORFLeft + nucleo #add new nucleus to ORF
        i+=3
    if ORFLeft != '':
        return [ORFLeft,i] #return orf and length
    else:
        return [None, i] #return None and length (solely to fit return style)


def find_all_ORFs_oneframe(dna):
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
    ORFs = []
    i=0
    while i<len(dna):
        newOrf = rest_of_ORF(dna[i:])
        if newOrf[0] == None: #if none then something went wrong or orf is done, end procees and return list of orfs
            return ORFs
        i+=newOrf[1]+3 #Skip over the orf AND THE STOP CODON
        ORFs.append(newOrf[0]) #add new orf to list of orfs
    return ORFs



def find_all_ORFs(dna):
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
    ORFs = []
    i = 0
    while i < 3:
        ORFs = ORFs + find_all_ORFs_oneframe(dna[i:]) #perform find orfs
        i+=1 #offset everything by 1 in order to then find nested orfs.
    return ORFs



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    otherRevDna = get_reverse_complement(dna)
    dnaORF = find_all_ORFs(dna)
    otherRevDnaOrf = find_all_ORFs(otherRevDna) #Do find all orfs for both forward dna and reverse compliment
    return dnaORF + otherRevDnaOrf #return list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    Orfs = find_all_ORFs_both_strands(dna)
    if not Orfs:
        return None
    lengs = []
    for Orf in Orfs:
        lengs.append(len(Orf)) #make list of orf lengs
    maxlen = max(lengs) #Find value of longest orf
    ind = lengs.index(maxlen) #Find index of longest ORF  IMPORTANT NOTE, ONLY RETURNS THE FIRST MAX FOUND, NOT ALL MAXS"""
    return Orfs[ind] #Return orf at said index


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i = 1
    ORFs = [longest_ORF(dna)]
    while i < num_trials:
        cutlocation = int(len(dna) * random.random()) #Pick random value within length of dna
        cutpiece = dna[:cutlocation] #split into piece 1
        restpiece = dna[cutlocation:] #and piece 2
        dna = restpiece + cutpiece #reverse order of the two pieces to choose another start location.
        ORFs.append(longest_ORF(dna)) #Append longest DNA of current orientation
        i+=1
    lengs = []
    for Orf in ORFs: #for each orf
        if Orf == None: #If something went wrong
            continue #skip this
        lengs.append(len(Orf)) #append lengths
    maxlen = max(lengs)
    ind = lengs.index(maxlen) #IMPORTANT NOTE, ONLY RETURNS THE FIRST MAX FOUND, NOT ALL MAXS"""
    return ORFs[ind]


def coding_strand_to_AA(dna, num_trials):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA", 10)
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT", 10)
        'MPA'
    """
    i = 0
    aaSequence = ''
    chain = longest_ORF_noncoding(dna, num_trials)
    while i < len(chain)-2:
        amino = aa_table[chain[i:i+3]]
        aaSequence += amino
        i+=3
    return aaSequence


def gene_finder(dna, num_trials):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    return coding_strand_to_AA(dna, num_trials)

def find_substring(substring, string):
    """
    >>> find_substring('ABABA', 'EIGHOSHONAGWABFNOGABABFNOABABANONGONOE')
    True
    >>> find_substring('WONDERFUL', 'WHAT A TERRIBLE WORLD')
    False
    """
    i = 0
    while i < len(string):
        if string[i] == substring[0]:
            j = 0
            k = 0
            while j < len(substring) and i + k < len(string):
                if string[i+k] != substring[j]:
                    break
                j+=1
                k+=1
            if j == len(substring):
                return True
        i+=1
    return False

def longest_common_substring(string1, string2, sizeAssume):
    """
    >>> longest_common_substring("What a wonderful world", "My life is a wonderful happening", 5)
    ' a wonderful '
    """
    results = []
    i=0
    while i < len(string1):
        j = 0
        if find_substring(string1[i:i+sizeAssume], string2):
            j = 1
            while j+i < len(string1):
                if find_substring(string1[i:i+sizeAssume+j], string2) == False or j+i == len(string1) - 1:
                    results.append(string1[i:i+sizeAssume+j-1])
                    break
                j+=1
        i+=j+1
    i=0
    lengs = []
    for res in results:
        lengs.append(len(res))
    maxlen = max(lengs)
    ind = lengs.index(maxlen)
    return results[ind]

def nitrogenase_finder(nitro, test_list, sizeAssume):
    '''
    >>> nitrogenase_finder('banana', [('test1', 'Holy shit, a banana'), ('test2', 'what a banan yo'), ('test3', 'I just got banned from the game')], 2)
    [('test1', 'banana'), ('test2', 'banan'), ('test3', 'ban')]
    '''
    results = []
    nitlen = len(nitro)
    i = 0
    testlen = len(test_list)
    for test in test_list:
        print(i + ' ' + testlen)
        if nitlen < len(test[1]):
            long_com_str = longest_common_substring(nitro, test[1], sizeAssume)
        else:
            long_com_str = longest_common_substring(test[1], nitro, sizeAssume)
        results.append((test[0], long_com_str))
    print(results)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    #dna = load_seq("./data/X73525.fa")
    #print(gene_finder(dna, 10000))

#print('This DNA is: EscU/YscU/HrcU family type III secretion system export apparatus switch protein [Salmonella enterica]')
