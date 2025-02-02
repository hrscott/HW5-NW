# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Define gap opening and gap extension penalties and substitution matrix file
    gap_open = -10
    gap_extend = -1
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"

    # Align all species to humans and print species in order of most similar to human BRD
    nw_hs = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    hs_scores = {}
    for seq, header in [(gg_seq, gg_header), (mm_seq, mm_header), (br_seq, br_header), (tt_seq, tt_header)]:
        _, _, alignB = nw_hs.align(hs_seq, seq)
        hs_scores[header] = nw_hs.alignment_score / len(hs_seq) # normalize by sequence length
    sorted_scores = sorted(hs_scores.items(), key=lambda x: x[1], reverse=True)
    print("Alignment scores relative to Homo sapiens BRD2:")
    for header, score in sorted_scores:
        print(header, score)

    # Print all of the alignment scores between each species BRD2 and human BRD2
    nw_gg = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    nw_mm = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    nw_br = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    nw_tt = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    hs_seqs = [hs_seq] * 4
    for seq, header, nw in [(gg_seq, gg_header, nw_gg), (mm_seq, mm_header, nw_mm), (br_seq, br_header, nw_br), (tt_seq, tt_header, nw_tt)]:
        _, _, alignB = nw.align(hs_seq, seq)
        print("Alignment score for {} and Homo sapiens BRD2: {}".format(header, nw.alignment_score))

if __name__ == "__main__":
    main()
