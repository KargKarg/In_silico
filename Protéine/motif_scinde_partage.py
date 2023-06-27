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
    return matrice


def initialiser_matrice(matrice: numpy.ndarray, sequence1: str, sequence2: str) -> numpy.ndarray:
    """
                    Permet d'obtenir la matrice de séquence commune.
                    Si sequence1[i-1] = sequence2[j-1]:
                        M[i,j] = M[i-1, j-1] + 1
                    Sinon:
                        M[i,j] = max(M[i, j-1], M[i-1, j])

                    Args:
                        - matrice (numpy.ndarray): Matrice remplie de 0.
                        - sequence1 (string): Première chaîne.
                        - sequence2 (string): Deuxième chaîne.

                    Returns:
                        - matrice (numpy.ndarray):  La matrice remplie pour la sous séquence.

    """
    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            if sequence1[i] == sequence2[j]:
                matrice[i+1, j+1] = matrice[i, j] + 1
            else:
                matrice[i+1, j+1] = max(matrice[i+1, j], matrice[i, j+1])
    return matrice


def motif_scinde_partage(fasta: str, enregistrer: bool = False) -> None:
    """
                Permet la lecture de la matrice de sous-séquence, et donne la sous-séquence scindée la plus longue.

                Args:
                    - fasta (string): Le chemin vers le fichier fasta.
                    - enregistrer (bool): Si sauvegarder = True enregistre dans un fichier texte.

                Returns:
                    - Aucun.

                Print:
                    - La sous-séquence commune la plus longue.

    """
    sequences = lire_fasta.lire(fasta)
    sequence1, sequence2 = sequences.values()
    matrice = initialiser_matrice(creer_matrice(sequence1, sequence2), sequence1, sequence2)
    texte = ''
    sous_seq = ''

    i, j = len(sequence1), len(sequence2)
    while i > 0 and j > 0:
        if sequence1[i-1] == sequence2[j-1]:
            sous_seq += sequence1[i-1]
            i, j = i-1, j-1
        elif matrice[i-1, j] > matrice[i, j-1]:
            i -= 1
        else:
            j -= 1

    sous_seq = sous_seq[::-1]
    texte += f"{list(sequences.keys())[0]} et {list(sequences.keys())[1]} partage ce motif réparti:\n\n{sous_seq}"

    if enregistrer:
        enregistrer_donnees.enregistrer("Résultats/motif_scinde_partage.txt", texte)

    print(texte)


motif_scinde_partage('Données/test.txt', False)
