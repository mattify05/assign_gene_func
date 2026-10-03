from Bio.Align import substitution_matrices
from Bio import SeqIO

cov_records = SeqIO.read("data/sars_cov_2.fa", "fasta")

print(len(cov_records.seq))

blosum62 = substitution_matrices.load("BLOSUM62")

d = 8

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)

    dp_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    dp_matrix[0][0] = 0

    # initialise all of first row
    for i in range(1, n + 1):
        dp_matrix[i][0] = dp_matrix[i - 1][0] + scoring_function(seq1[i-1], "-")

    for j in range(1, m + 1):
        dp_matrix[0][j] = dp_matrix[0][j - 1] + scoring_function("-", seq2[j - 1])

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = dp_matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = dp_matrix[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = dp_matrix[i][j - 1] + scoring_function("-", seq2[j - 1])

            dp_matrix[i][j] = max(diagonal, up, left)

    i = n
    j = m

    aligned_seq1 = []
    aligned_seq2 = []

    while i > 0 or j > 0:

        if i > 0 and j > 0:
            diagonal = dp_matrix[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])

            if dp_matrix[i][j] == diagonal:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])

                i -= 1
                j -= 1
                continue

        if i > 0:

            up = dp_matrix[i - 1][j] + scoring_function(seq1[i - 1], "-")

            if dp_matrix[i][j] == up:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append("-")

                i -= 1

                continue

        if j > 0:

            left = dp_matrix[i][j - 1] + scoring_function("-", seq2[j - 1])

            if dp_matrix[i][j] == left:
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[j - 1])

                j -= 1
                continue

    aligned_sequence1 = "".join(aligned_seq1[::-1])
    aligned_sequence2 = "".join(aligned_seq2[::-1])

    return aligned_sequence1, aligned_sequence2, float(dp_matrix[n][m])


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    return 0

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == "-" or aa_j == "-":
        return -d
    
    return float(blosum62[aa_i, aa_j])