from Module import lire_fasta, enregistrer_donnees


def epissage(fasta: str, sauvegarder: bool = False) -> None:
    """
            Prend un fichier FASTA, avec une séquence au début et une suite de séquences qui seront des introns.
            Supprime sur la séquence d'intérêt la succession d'intron.

            Args:
                - fasta (str): La chaîne de caractères qui désigne le chemin du fichier.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                Aucun.

            Print:
                - La séquence après épissage.

    """
    seq = lire_fasta.lire(fasta)
    introns = iter(seq.values())
    sequence = introns.__next__()
    for intron in introns:
        sequence = sequence.split(intron)[0] + sequence.split(intron)[1]

    texte = f"{iter(seq.keys()).__next__()}\n{sequence}"
    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/epissage_résultats.txt', texte)

    print(sequence)


epissage('Données/epissage.txt', True)
