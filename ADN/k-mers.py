from Module import table_kmers, lire_fasta, enregistrer_donnees


def kmers(fasta: str, k: int = 4, sauvegarder: bool = False):
    """
        Prend un chemin vers le fichier FASTA, un entier K correspondant à la taille du k-mers.
        Affiche pour chaque séquence de ce fichier FASTA, le nombre de combinaison du k-mers retrouvé.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - k (int): L'entier qui correspond à la fenêtre du k-mers.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Le résultat de la recherche de k-mers.

    """
    seqs = lire_fasta.lire(fasta)
    tableau = table_kmers.kmers(k)

    texte = ''

    for id, sequence in seqs.items():
        texte += f'{id}\n'
        for i in range(len(sequence)-k+1):
            tableau[sequence[i:i+k]] += 1
        texte += '\n'.join([f"{key} {val}" for key, val in tableau.items()])
        tableau = table_kmers.kmers(k)

    if sauvegarder:
        enregistrer_donnees.enregistrer(f'Résultats/{k}mers_résultats.txt', texte)

    print(texte)

kmers('Données/test.txt', 4, True)