from In_silico.Module import lire_fasta, enregistrer_donnees
import numpy as np


def p_distance(sequence1, sequence2):
    difference = 0
    for i in range(min(len(sequence1), len(sequence2))):
            if sequence1[i] != sequence2[i]:
                difference += 1
    return round(difference/min(len(sequence1), len(sequence2)), 5)
def matrice_distance(fasta, sauvegarder):
    sequences = lire_fasta.lire(fasta)
    taille = len(sequences.values())
    matrice = np.full((taille, taille), 0, float)
    texte = ""
    for i in range(len(sequences.values())):
        for j in range(len(sequences.values())):
            matrice[i, j] += p_distance(list(sequences.values())[i], list(sequences.values())[j])
            texte += f"{round(matrice[i,j], 5)} "
        texte += '\n'
    print(texte)

matrice_distance('Données/test.txt', 1)