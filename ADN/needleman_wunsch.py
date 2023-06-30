from In_silico.Module import lire_fasta, enregistrer_donnees
import numpy as np


def initialisation_matrice(sequence1: str, sequence2: str, identite: int, substitution: int, gap: int) -> tuple:
    """
            Permet de créer une matrice en fonction de la taile des séquences, avec séquence1 de taille n,
            et séquence2 de taille m.
            La matrice sera de taille (n,m).
            Et applique l'algorithme de Needleman et Wunsch pour la remplir.

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

    for i in range(1, lignes):
        matrice[i, 0] = matrice[i-1, 0] + gap

    for i in range(1, colonnes):
        matrice[0, i] = matrice[0, i-1] + gap

    for i in range(1, lignes):
        for j in range(1, colonnes):
            haut_gauche = matrice[i-1, j-1]
            gauche = matrice[i, j-1]
            haut = matrice[i-1, j]
            if sequence1[i-1] == sequence2[j-1]:
                matrice[i, j] = max(haut_gauche+identite, gauche+gap, haut+gap)
            else:
                matrice[i, j] = max(haut_gauche+substitution, gauche+gap, haut+gap)

    return lecture(matrice, sequence1, sequence2, identite, substitution, gap)


def lecture(matrice: np.ndarray, sequence1: str, sequence2: str, identite: int, substitution: int, gap: int) -> tuple:
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
    seq1 = ''
    evenement = ''
    seq2 = ''

    cpt_identite = 0
    cpt_substitution = 0
    cpt_gap = 0

    i, j = len(sequence1), len(sequence2)

    while i > 0 and j > 0:

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
                evenement += ' '
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

    texte += seq1[::-1] + '\n' + evenement[::-1] + '\n' + seq2[::-1]

    return matrice, texte, matrice[len(sequence1), len(sequence2)],\
        cpt_identite, cpt_substitution, cpt_gap, max(len(seq1), len(seq2))


def needleman_wunsch(fasta, identite, substitution, gap, sauvergarder):
    """
            Prend un fasta avec 2 séquences, et applique à travers d'autres fonctions l'algorithme de
            Needleman et Wunsch avec les valeurs données par l'utilisateur.

            Args:
                - fasta (string): Le chemin vers le fichier fasta.
                - identite (int): Le score correspondant à un événement d'identité.
                - substitution (int): Le score correspondant à un événement de substitution.
                - gap (int): Le score correspondant à un événement d'indel.
                - enregistrer (bool): Si sauvegarder = True enregistre dans un fichier texte.

            Returns:
                - Aucun.

            Print:
                - L'alignement des deux séquences, les métriques de l'alignement et la matrice associée.

    """
    sequences = lire_fasta.lire(fasta)
    sequence1, sequence2 = sequences.values()
    sid1, sid2 = sequences.keys()
    texte = ''
    mat = ""

    matrice, align, score, cid, csub, cgap, taille = \
        initialisation_matrice(sequence1, sequence2, identite, substitution, gap)

    texte += f"{sid1}\n{align}\n{sid2}\n\n\nScore: {score}\n"
    texte += f"Identite: {cid}/{taille} ({round(cid/taille*100, 2)}%)\n"
    texte += f"Substitution: {csub}/{taille} ({round(csub / taille * 100, 2)}%)\n"
    texte += f"Gap: {cgap}/{taille} ({round(cgap / taille * 100, 2)}%)\n\n"
    texte += f"Valeurs:\nIdentite: {identite}\nSubstitution: {substitution}\nGap: {gap}"

    for i in range(len(sequence1)+1):
        for j in range(len(sequence2)+1):
            mat += f"{matrice[i, j]} "
        mat += '\n'

    if sauvergarder:
        enregistrer_donnees.enregistrer('Résultats/matrice_needleman_wunsch.txt', mat)
        enregistrer_donnees.enregistrer('Résultats/score_needleman_wunsch.txt', texte)

    print(texte)


needleman_wunsch('Données/test.txt', 2, -1, -2, True)
