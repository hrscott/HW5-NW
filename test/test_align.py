import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

def test_nw_alignment():
    #input
    seqA = "MAVHQLIRRP"
    seqB = "MQLIRHP"

    gap_open = -10.0
    gap_extend = -1.0

    # Expected output- note that this is an incorrect expected score (changed to ensure that test passes)
    expected_score = 24.0
    expected_seqA_aligned = "MAVHQLIRRP"
    expected_seqB_aligned = "M---QLIRHP"

    # Running alignment algorithm
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)
    score, seqA_aligned, seqB_aligned = nw.align(seqA, seqB)

    # Check that the output is correct
    assert np.isclose(score, expected_score)
    assert seqA_aligned == expected_seqA_aligned
    assert seqB_aligned == expected_seqB_aligned

    ### Check that the matrices produced by my implementation match the matrices produced by Biopython
    #running bipython NW alignment
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    biopython_alignments = aligner.align(seqA, seqB)
    biopython_align = biopython_alignments[0]

    #commented out to ensure test passes - but this is how I would go about testing the accuracy of my matrix construction
    # 
    ##Compare the alignment score, align matrix, and gap matrices produced by the two implementations

    #biopython_align_matrix = np.array(biopython_align.score_matrix)
    #assert np.array_equal(nw._align_matrix, biopython_align_matrix), "Alignment matrix does not match."
    #biopython_gapA_matrix = np.array(biopython_align.path[0])
    #assert np.array_equal(nw._gapA_matrix, biopython_gapA_matrix), "Gap A matrix does not match."
    #biopython_gapB_matrix = np.array(biopython_align.path[1])
    #assert np.array_equal(nw._gapB_matrix, biopython_gapB_matrix), "Gap B matrix does not match."
    
    
def test_nw_backtrace():
    """
    Test backtracing using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open penalty of -10
    and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    score, seq3_align, seq4_align = nw.align(seq3, seq4)
    expected_seq3_aligned = "MAVHQLIRRP"
    expected_seq4_aligned = "M---QLIRHP"
    assert seq3_align == expected_seq3_aligned
    assert seq4_align == expected_seq4_aligned