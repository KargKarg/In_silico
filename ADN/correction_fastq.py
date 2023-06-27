from In_silico.Module import lire_fastq, enregistrer_donnees, table_phred


def moyenne(sequence: str) -> float:
    """
                Prend une séquence qualité et renvoie la moyenne de la qualité.

                Args:
                    - sequence (str): La chaîne de caractères correspondant à la séquence qualité.

                Returns:
                    - moy (float): La moyenne qualité de la séquence.

                Print:
                    - Aucun.

    """
    moy = 0
    scores = table_phred.table()
    for symbole in sequence:
        moy += scores[symbole]
    moy /= len(sequence)
    return moy


def eroder(sequence: str, seuil: int, mini: int) -> tuple:
    """
                Prend une séquence qualité et calcul la moyenne en érodant l'extrémité finale.
                Renvoie un tuple qui sera un indice et la moyenne de la dernière séquence.
                Cela se base sur le postulat que les séquençages type ""Illumina"" ont une qualité décroissante.
                L'algorithme ne descendra pas en dessous de la valeur imposée.

                Args:
                    - sequence (str): La chaîne de caractères correspondant à la séquence qualité.
                    - seuil (int): L'entier qui correspond à la limite de valabilité.
                    - mini (int): La taille minimale requise pour qu'une séquence soit considérée.

                Returns:
                    - stop (int): Correpond à l'indice où la moyenne >= seuil.
                    - (float): La qualité moyenne de la séquence.

                Print:
                    - Aucun.

    """
    stop = 0
    for i in range(len(sequence), mini-2, -1):
        stop = i
        if moyenne(sequence[:i]) >= seuil:
            break
    return stop, moyenne(sequence[:stop])


def correction(fastq, seuil, mini, sauvegarder) -> None:
    """
                Prend un chemin vers le FASTQ, un seuil qui déterminera la valeur considérée et une valeur minimum
                requise pour qu'une séquence soit aussi considérée.
                Permet de discriminer l'ensemble des séquences en fonction de leur qualité moyenne,
                 en vérifiant si elles dépassent ou non le seuil établi.

                Args:
                    - fastq (str): La chaîne de caractères du chemin vers le fichier.
                    - seuil (int): L'entier qui correspond à la limite de viabilité.
                    - mini (int): L'entier qui donne la taille limite de viabilité.
                    - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

                Returns:
                    Aucun.

                Print:
                    - Les séquences sont enregistrés dans Données et le log dans Résultats.

    """
    sequences = lire_fastq.lire(fastq)
    nb_reads = len(sequences.values())
    mauvaises = []
    bonnes_sequences = ''
    mauvaises_sequences = ''
    texte_supp = ''
    texte_inf = ''
    cpt = 0

    for sid, sequence in sequences.items():

        if moyenne(sequence[1]) < seuil:

            indice, moy = eroder(sequence[1], seuil, mini)

            if moy >= seuil:
                texte_supp += f"{sid}:\nQualite: {round(moyenne(sequence[1]), 2)}\nTaille: {len(sequence[1])}pb\n\n"
                texte_supp += f"Apres erodage:\nQualite: {round(moy, 2)}\nTaille: {indice}pb\n\n\n\n"
                sequences[sid][0], sequences[sid][1] = sequences[sid][0][:indice], sequences[sid][1][:indice]

            else:
                cpt += 1
                texte_inf += f"{sid}:\nQualite: {round(moyenne(sequence[1]), 2)}\nTaille: {len(sequence[1])}pb\n\n"
                texte_inf += f"Apres erodage:\nQualite: {round(moy, 2)}\nTaille: {indice}pb\nElle n'est pas utilisable car inferieur a la taille minimum\n\n\n\n"
                mauvaises.append(sid)

    for sid in mauvaises:
        mauvaises_sequences += f"{sid}\n{sequences[sid][0]}\n+\n{sequences[sid][1]}\n"
        del sequences[sid]

    for sid, sequence in sequences.items():
        bonnes_sequences += f"{sid}\n{sequence[0]}\n+\n{sequence[1]}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/logs_erodage.txt', f"{cpt} reads non conformes\n{nb_reads-cpt} reads conformes\n\n{texte_inf}")
        enregistrer_donnees.enregistrer('Données/bonnes_sequences_erodage.txt', bonnes_sequences)
        enregistrer_donnees.enregistrer('Données/mauvaises_sequences_erodage.txt', mauvaises_sequences)


correction('Données/reads1_fastq.txt', 20, 25, True)
