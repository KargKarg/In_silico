from Module import lire_fasta, enregistrer_donnees


def super_string(fasta: str, sauvegarder: bool = False):
    seqs = lire_fasta.lire(fasta)


super_string("Données/test.txt")
