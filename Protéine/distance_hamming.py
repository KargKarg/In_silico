from In_silico.Module import lire_fasta, enregistrer_donnees


def distance_hamming(fasta: str, sauvegarder: bool = False) -> None:
    """
                Prend un chemin vers le fichier FASTA contenant 2 séquences.
                Et sauvegarder qui permet de demander à l'utilisateur si il veut enregistrer le résulat.
                Permet d'obtenir la distance Hamming entre 2 séquences en la parcourant tout le long.


                Args:
                    - fasta (str): La chaîne de caractères du chemin vers le fichier.
                    - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

                Returns:
                    - Aucun.

                Print:
                    - La distance Hamming avec les points de mutations.

    """
    sequence1, sequence2 = lire_fasta.lire(fasta).values()
    ids1, ids2 = lire_fasta.lire(fasta).keys()
    distance = 0
    indices = []
    texte = f"{ids1}\n"

    for i in range(min(len(sequence1), len(sequence2))):
        if sequence1[i] != sequence2[i]:
            distance += 1
            indices.append(i)

    texte += f'{sequence1}\n'

    for i in range(min(len(sequence1), len(sequence2))):
        if i in indices:
            texte += 'x'
        else:
            texte += ' '
    texte += f'\n{sequence2}\n{ids2}\n\nDistance Hamming de {distance}'

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/distance_hamming.txt', texte)

    print(texte)


distance_hamming('Données/test.txt', True)
