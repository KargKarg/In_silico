from Module import lire_fasta, enregistrer_donnees


def proportion(fasta: str, sauvegarder: bool = False) -> None:
    """
        Prend un chemin vers le fichier FASTA.
        Affiche pour chaque séquence du FASTA, les proportions en AT et GC.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - La proportion en AT et GC.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ""
    base_comp = ('AT', 'GC')
    for id, sequence in seqs.items():
        data = [(bases, (sequence.count(bases[0])+sequence.count(bases[1]))/(len(sequence))*100) for bases in base_comp]
        texte += f"{id}\n{data}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/proportion_résultats.txt', texte)

    print(texte)


proportion("Données/test.txt", True)
