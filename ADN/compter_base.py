from Module import lire_fasta, enregistrer_donnees


def compter_base(fasta: str, sauvegarder: bool = False):
    """
        Prend un chemin vers le fichier FASTA.
        Affiche pour chaque séquence de ce fichier FASTA, le nombre de chaque nucléotide qu'il posséde.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Le résultat du compte des bases.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ''
    for id, sequence in seqs.items():
        data = [(base, sequence.count(base)) for base in set(sequence)]
        texte += f"{id}\n{data}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/compte_base_résultats.txt', texte)

    print(texte)


compter_base('Données/seq.txt')
