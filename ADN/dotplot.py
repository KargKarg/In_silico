from In_silico.Module import lire_fasta


def dotplot(fasta: str, fenetre: int = 3, sauvegarder: bool = False):
    sequences = lire_fasta.lire(fasta)
    seq1, seq2 = sequences.values()[0], sequences.values()[1]

    points = [0 for _ in range(max(len(seq1), len(seq2)))]

    for i in range(0, len(points)):
        if seq1[i:i+fenetre] == seq2[i:i+fenetre]:
            points[i] += 1
