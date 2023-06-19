from Module import masse_aa, lire_fasta, enregistrer_donnees


def masse_prot(fasta: str, sauvegarder: bool = False):
    """
            Prend un fichier FASTA.
            Affiche pour chaque séquence du FASTA, sa masse monoisotopic.

            Args:
                - fasta (str): La chaîne de caractères qui désigne le chemin du fichier.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                Aucun.

            Print:
                - La masse des protéines.

    """
    tableau = masse_aa.tableau_masse()
    seqs = lire_fasta.lire(fasta)
    texte = ""
    for id, sequence in seqs.items():
        poids = 0
        for AA in sequence:
            poids += tableau[AA]
        texte += f"{id}\n{poids}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/masse_proteine_résultats.txt', texte)

    print(texte)


masse_prot('Données/test.txt', True)
