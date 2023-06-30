import numpy
import numpy as np
from In_silico.Module import lire_fasta, enregistrer_donnees


def creer_matrice(sequence1: str, sequence2: str) -> numpy.ndarray:
    """
                        Permet de créer une matrice en fonction de la taile des séquences, avec séquence1 de taille n,
                        et séquence2 de taille m.
                        La matrice sera de taille (n,m).

                        Args:
                            - sequence1 (string): Première chaîne.
                            - sequence2 (string): Deuxième chaîne.

                        Returns:
                            - matrice (numpy.ndarray):  La matrice.

    """
    matrice = np.full((len(sequence1)+1, len(sequence2)+1), 0)

    for i in range(len(sequence1)+1):
        matrice[i, 0] = i

    for j in range(len(sequence2)+1):
        matrice[0, j] = j

    return matrice


def distance_de_modification(fasta: str, sauvergarder: bool = False) -> None:
    """
                    Permet d'avoir la distance de modification entre 2 séquences, c'est à dire le nombre
                    d'indel ou de substitution pour transformer la séquence1 en séquence2.

                    Args:
                        - fasta (string): Le chemin vers le fichier fasta.
                        - enregistrer (bool): Si sauvegarder = True enregistre dans un fichier texte.

                    Returns:
                        - Aucun.

                    Print:
                        - La matrice avec la distance de modification.

    """
    sequences = lire_fasta.lire(fasta)
    sequence1, sequence2 = sequences.values()
    matrice = creer_matrice(sequence1, sequence2)
    texte = ""
    mat = ""

    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            gauche = matrice[i-1, j] + 1
            bas = matrice[i, j-1] + 1
            gauche_bas = matrice[i-1, j-1]
            if sequence1[i - 1] != sequence2[j - 1]:
                gauche_bas += 1
            matrice[i, j] = min(gauche, bas, gauche_bas)

    for i in range(len(sequence1)+1):
        for j in range(len(sequence2)+1):
            mat += f"{matrice[i, j]} "
        mat += '\n'

    texte += f"La distance de modification est de {matrice[len(sequence1), len(sequence2)]}"

    if sauvergarder:
        enregistrer_donnees.enregistrer('Résultats/modification_distance.txt', texte)
        enregistrer_donnees.enregistrer('Résultats/matrice_modification_distance.txt', mat)

    print(texte)


distance_de_modification("Données/test.txt", True)
