from Module import lire_fasta, enregistrer_donnees


def graphe(fasta: str, nb: int = 5, sauvegarder: bool = False) -> None:
    """
        Prend un chemin vers le FASTA fasta, un entier nb correspondant à la taille du chevauchement considéré.
        Affiche pour chaque séquence de ce ficher FASTA, les arcs (u,v) où :

                                                - u est le fragment finissant par le motif.
                                                - v est le fragment commençant par le motif.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - nb (int): L'entier qui correspond à la taille du chevauchement considéré.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Les arcs du graphe.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ""
    for Fid, Fsequence in seqs.items():
        for Sid, Ssequence in seqs.items():
            if Fid != Sid and Fsequence[-nb:] == Ssequence[:nb]:
                texte += f"{Fid[1:]}--->{Sid[1:]}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/graphe_chevauchant_résultats.txt', texte)

    print(texte)


graphe('Données/test.txt', 3, True)
