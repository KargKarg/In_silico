from In_silico.Module import lire_fasta


def super_string(fasta: str, sauvegarder: bool = False):
    seqs = lire_fasta.lire(fasta)


super_string("Données/test.txt")
