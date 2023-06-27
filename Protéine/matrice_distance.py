from In_silico.Module import lire_fasta, enregistrer_donnees, p_distance
import numpy as np


def matrice_distance(fasta: str, sauvegarder: bool = False) -> None:
    """
                Permet d'obtenir la matrice des ""p distances"" entre chacunes des séquences.
                M[i,j] = ""p distances"" entre i et j.

                Args:
                    - chaine1 (string): Première chaîne.
                    - chaine2 (string): Deuxième chaîne.

                Returns:
                    - matrice (numpy.ndarray):  La matrice des ""p distances"".

    """
    sequences = lire_fasta.lire(fasta)
    taille = len(sequences.values())
    matrice = np.full((taille, taille), 0, float)
    texte = ""
    indices = "Pour les indices on a :\n"

    for i in range(len(sequences.values())):
        indices += f"{i+1} --> {list(sequences.keys())[i]}\n"

        for j in range(len(sequences.values())):

            matrice[i, j] += p_distance.calcul(list(sequences.values())[i], list(sequences.values())[j])
            texte += f"{round(matrice[i,j], 5)} "

        texte += '\n'

    texte += f"\n{indices}"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/matrice_distance.txt', texte)

    print(texte)
    return matrice


matrice_distance('Données/test.txt', True)
