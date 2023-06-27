from In_silico.Module import lire_fasta, enregistrer_donnees
import numpy as np


def profile(fasta: str, sauvegarder: bool = False) -> None:
    """
        Prend un chemin vers le fichier FASTA.
        Affiche la séquence profile avec l'aide d'une matrice.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - La séquence profile.

    """
    seqs = lire_fasta.lire(fasta)
    matrice = np.empty((len(seqs), len(max(seqs.values(), key=lambda x: len(x[1])))), dtype=object)
    cpt = 0
    for sequence in seqs.values():
        for i in range(len(sequence)):
            matrice[cpt, i] = sequence[i]
        cpt += 1
    sequence_profile = ''
    for i in range(len(matrice.T)):
        maximum = 0
        nt = ''
        for elem in matrice.T[i]:
            if elem is not None and matrice.T[i].tolist().count(elem) > maximum:
                maximum = matrice.T[i].tolist().count(elem)
                nt = elem
        sequence_profile += nt

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/profile_sequence_résultats.txt', sequence_profile)

    print(sequence_profile)


profile('Données/test.txt', True)
