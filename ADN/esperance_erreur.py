from In_silico.Module import fastq_erreur, lire_fastq, enregistrer_donnees


def esperance(fastq, sauvegarder):
    """
                Prend un chemin vers le fichier FASTQ.
                Calcul l'esperance d'erreur sur chaque séquence dans le FASTQ.

                Args:
                    - fastq (str): La chaîne de caractères du chemin vers le fichier.

                Returns:
                    - Aucun.

                Print:
                    - L'espérance pour chaque séquence dans le FASTQ.

    """
    sequences = lire_fastq.lire(fastq)
    texte = ''

    for sid, sequence in sequences.items():
        valeur = 0

        for symbole in sequence[1]:
            valeur += fastq_erreur.proba(symbole)

        texte += str(sid) + ' ' + str(round(valeur, 2)) + ' ' + str(len(sequence[0])) + '\n'

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/esperance_erreur_sequences.txt', texte)

    print(texte)


esperance('Données/mauvaises_sequences_erodage.txt', True)
