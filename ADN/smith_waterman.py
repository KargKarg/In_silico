import numpy as np
from Module import lire_fasta, enregistrer_donnees


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
                - matrice (numpy.ndarray): Matrice remplie par l'algorithme de Smith et Waterman.
                - maximum[1] (tuple): Correspond à la position de la valeur maximum dans la matrice.

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

    return matrice, maximum[1]


def lecture(matrice: np.ndarray, lignes: int, colonnes: int, sequence1: str, sid1: str, sequence2: str, sid2: str, identite: int, substitution: int, gap: int) -> str:
    """
            Permet d'avoir le score d'alignement entre 2 séquences, ainsi que l'alignement, le nombre d'identite,
            de substitution et de gap.
            Cette fonction prend une matrice déjà complète.

            Args:
                - matrice (numpy.ndarray): La matrice complète.
                - lignes (int): Correspond à la ligne de la valeur maximum de la matrice.
                - colonnes (int): Correspond à la colonne de la valeur maximum de la matrice.
                - sequence1 (string): Première chaîne.
                - sid1 (string): Correspond à l'id de la séquence 1.
                - sequence2 (string): Deuxième chaîne.
                - sid2 (string): Correspond à l'id de la séquence 2.
                - identite (int): Le score correspondant à un événement d'identité.
                - substitution (int): Le score correspondant à un événement de substitution.
                - gap (int): Le score correspondant à un événement d'indel.

            Returns:
                - texte (string): L'alignement ainsi que les métriques formatés.

    """
    texte = ''
    evenement = ''
    seq1 = f"{sequence1[lignes:][::-1]}"
    seq2 = f"{sequence2[colonnes:][::-1]}"

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

    seq1 += f"{sequence1[:i][::-1]}"
    seq2 += f"{sequence2[:j][::-1]}"

    for _ in range(max(lignes, colonnes)-len(evenement)):
        evenement += " "

    for _ in range(max(lignes, colonnes)-len(seq1)):
        seq1 += " "

    for _ in range(max(lignes, colonnes)-len(seq2)):
        seq2 += " "

    texte += f"{sid1}\n{seq1[::-1]}\n{evenement[::-1]}\n{seq2[::-1]}\n{sid2}\n\n\n"

    texte += f"Identite: {cpt_identite}/{min(len(seq1.replace(' ', '')), len(seq2.replace(' ', '')))} ({round(cpt_identite / min(len(seq1.replace(' ', '')), len(seq2.replace(' ', ''))) * 100)}%)\n"
    texte += f"Substitution: {cpt_substitution}/{min(len(seq1.replace(' ', '')), len(seq2.replace(' ', '')))} ({round(cpt_substitution / min(len(seq1.replace(' ', '')), len(seq2.replace(' ', ''))) * 100)}%)\n"
    texte += f"Gap: {cpt_gap}/{min(len(seq1.replace(' ', '')), len(seq2.replace(' ', '')))} ({round(cpt_gap / min(len(seq1.replace(' ', '')), len(seq2.replace(' ', ''))) * 100)}%)\n\n\n"

    texte += f"Avec:\nIdentite: {identite}\nSubstitution: {substitution}\nGap: {gap}"

    return texte


def smith_waterman(fasta: str, identite: int, substitution: int, gap: int, sauvegarder: bool = False) -> None:
    """
        Prend un fasta avec 2 séquences, et applique à travers d'autres fonctions l'algorithme de
        Smith et Waterman avec les valeurs données par l'utilisateur.

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
    mat = ''

    matrice, coords = initialisation_matrice(sequence1, sequence2, identite, substitution, gap)
    x, y = coords
    texte = lecture(matrice, x, y, sequence1, sid1, sequence2, sid2, identite, substitution, gap)

    for i in range(len(sequence1)+1):
        for j in range(len(sequence2)+1):
            mat += f"{matrice[i, j]} "
        mat += '\n'

    if sauvegarder:
        enregistrer_donnees.enregistrer("Résultats/smith_waterman.txt", texte)
        enregistrer_donnees.enregistrer("Résultats/matrice_smith_waterman.txt", mat)

    print(texte)


smith_waterman("Données/test.txt", 2, 0, -1, True)

