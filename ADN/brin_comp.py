from In_silico.Module import lire_fasta, enregistrer_donnees


def brin_complementaire(fasta: str, sauvegarder: bool = False) -> dict:
    """
        Prend un chemin vers le fichier FASTA..
        Affiche pour chaque séquence de ce fichier FASTA, le brin complémentaire de 5' à 3'.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            - seqs_comp (dict): Un dictionnaire avec {id de la séquence: brin complémentaire de la séquence}.

        Print:
            - Le brin complémentaire.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ''
    seqs_comp = {}
    for id, sequence in seqs.items():
        sequence = sequence.replace('T', 'a').replace('A', 't').replace('G', 'c').replace('C', 'g')
        seqs_comp[id] = sequence[::-1].upper()
        texte += f"{id}\n{sequence[::-1].upper()}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/brin_complémentaire_résultats_5a3.txt', texte)

    return seqs_comp


brin_complementaire('Données/seq.txt')
