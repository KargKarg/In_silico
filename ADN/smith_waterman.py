import numpy as np


def initialisation_matrice(sequence1: str, sequence2: str, identite: int, substitution: int, gap: int) -> tuple:
    """
            Permet de créer une matrice en fonction de la taile des séquences, avec séquence1 de taille n,
            et séquence2 de taille m.
            La matrice sera de taille (n,m).
            Et applique l'algorithme de Smith et Waterman pour la remplir.

            Args:
                - sequence1 (string): Première chaîne.
                - sequence2 (string): Deuxième chaîne.
                - identite (int): Le score correspondant à un événement d'identité.
                - substitution (int): Le score correspondant à un événement de substitution.
                - gap (int): Le score correspondant à un événement d'indel.

            Returns:
                - lecutre (function): Fonction qui donne les métriques et l'alignement.

    """
    lignes, colonnes = len(sequence1)+1, len(sequence2)+1
    matrice = np.full((lignes, colonnes), 0)
    maximum = [0, ()]

    for i in range(1, lignes):
        for j in range(1, colonnes):
            haut_gauche = matrice[i-1, j-1]
            gauche = matrice[i, j-1]
            haut = matrice[i-1, j]
            if sequence1[i-1] == sequence2[j-1]:
                matrice[i, j] = max(haut_gauche+identite, gauche+gap, haut+gap, 0)
            else:
                matrice[i, j] = max(haut_gauche+substitution, gauche+gap, haut+gap, 0)
            if matrice[i, j] > maximum[0]:
                maximum[0] = matrice[i, j]
                maximum[1] = (i, j)

    lecture(matrice, maximum[1][0], maximum[1][1], sequence1, sequence2, identite, substitution, gap)

    return matrice, maximum


def lecture(matrice: np.ndarray, lignes: int, colonnes: int, sequence1: str, sequence2: str, identite: int, substitution: int, gap: int) -> tuple:
    """
            Permet d'avoir le score d'alignement entre 2 séquences, ainsi que l'alignement, le nombre d'identite,
            de substitution et de gap.
            Cette fonction prend une matrice déjà complète.

            Args:
                - matrice (numpy.ndarray): La matrice complète.
                - sequence1 (string): Première chaîne.
                - sequence2 (string): Deuxième chaîne.
                - identite (int): Le score correspondant à un événement d'identité.
                - substitution (int): Le score correspondant à un événement de substitution.
                - gap (int): Le score correspondant à un événement d'indel.

            Returns:
                - matrice (numpy.ndarray): La matrice complète.
                - texte (string): L'alignement ainsi que les métriques formatés.
                - : Le score de l'alignement
                - cpt_identite (int): Le nombre d'identité.
                - cpt_substitution (int): Le nombre d'identité.
                - cpt_gap (int): Le nombre d'identité.
                 - (int): La taille de l'alignement.

    """
    texte = ''
    evenement = ''
    seq1 = f"{sequence1[lignes:]}"
    seq2 = f"{sequence2[colonnes:]}"

    cpt_identite = 0
    cpt_substitution = 0
    cpt_gap = 0

    i, j = lignes, colonnes

    while i > 0 and j > 0 and matrice[i, j] > 0:

        if sequence1[i-1] == sequence2[j-1]:
            haut_gauche = matrice[i - 1, j - 1]+identite
            gauche = matrice[i, j - 1]+gap
            haut = matrice[i - 1, j]+gap
            direction = max(haut_gauche, gauche, haut)
        else:
            haut_gauche = matrice[i - 1, j - 1]+substitution
            gauche = matrice[i, j - 1]+gap
            haut = matrice[i - 1, j]+gap
            direction = max(haut_gauche, gauche, haut)

        if haut_gauche == direction:
            seq1 += sequence1[i-1]
            seq2 += sequence2[j-1]
            if sequence1[i-1] == sequence2[j-1]:
                evenement += '|'
                cpt_identite += 1
            else:
                evenement += ':'
                cpt_substitution += 1
            i, j = i-1, j-1

        elif haut == direction:
            seq1 += sequence1[i-1]
            seq2 += '-'
            evenement += ' '
            i -= 1
            cpt_gap += 1

        else:
            seq1 += '-'
            seq2 += sequence2[j-1]
            evenement += ' '
            j -= 1
            cpt_gap += 1

    seq1 += f"{sequence1[:i]}"
    seq2 += f"{sequence1[:j]}"

    print(seq1[::-1], '\n', evenement[::-1], '\n', seq2[::-1])

print(initialisation_matrice('GACTTAC', 'CGTGAATTCAT', 2, -1, -2))
