# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import operator

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

    nw = NeedlemanWunsch(sub_matrix_file='BLOSUM62', gap_open=-10, gap_extend=-1)

    # Align all species to humans and print species in order of most similar to human BRD
    species = {"Homo sapiens": hs_seq, "Gallus gallus": gg_seq, "Mus musculus": mm_seq, "Balaeniceps rex": br_seq, "Tursiops truncatus": tt_seq}
    species_alignment_scores = {}
    for sp_name, sp_seq in species.items():
        alignment_score, _, _ = nw.align(hs_seq, sp_seq)
        species_alignment_scores[sp_name] = alignment_score

    sorted_species = sorted(species_alignment_scores.items(), key=operator.itemgetter(1), reverse=True)
    print("Species in order of most similar to human BRD2:")
    for sp_name, sp_score in sorted_species:
        print(sp_name)

    # Print all of the alignment score between each species BRD2 and human BRD2
    print("\nAlignment scores between each species BRD2 and human BRD2:")
    for sp_name, sp_seq in species.items():
        alignment_score, _, _ = nw.align(hs_seq, sp_seq)
        print("{}: {}".format(sp_name, alignment_score))
    

if __name__ == "__main__":
    main()