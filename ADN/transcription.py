from In_silico.Module import lire_fasta, enregistrer_donnees


def transcription(fasta: str, sauvegarder: bool = False) -> dict:
    """
        Prend un chemin vers le fichier FASTA.
        Affiche pour chaque séquence du fichier FASTA, la transcription du brin.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            - sequences (dict): Un dictionnaire avec {id de la séquence: brin complémentaire de la séquence}.

        Print:
            - La transcription des brins.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ''
    sequences = {}
    for id, sequence in seqs.items():
        sequence = sequence[::-1].replace('A', 'u').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
        texte += f"{id}\n{sequence}\n"
        sequences[id[1:]] = sequence

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/transcription_résultats_5a3.txt', texte)

    return sequences


transcription('Données/test.txt', True)
