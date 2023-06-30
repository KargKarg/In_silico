from In_silico.Module import lire_fasta, enregistrer_donnees


def hamming(sequence1: str, sequence2: str) -> int:
    """
                    Donne la distance Hamming entre les deux séquences.

                    Args:
                        - sequence1 (str): La chaîne de caractères qui correspond à la séquence une.
                        - sequence2 (str): La chaîne de caractères qui correspond à la séquence deux.

                    Returns:
                        - hmg (int): La distance Hamming calculée.

                    Print:
                        - Aucun.

    """
    hmg = 0
    for i in range(min(len(sequence1), len(sequence2))):
        if sequence1[i] != sequence2[i]:
            hmg += 1
    return hmg


def pattern_approximatif(fasta: str, motif: str, mismatch: int, sauvegarder: bool = False) -> None:
    """
                    Prend un chemin vers le fichier FASTA contenant 2 séquences.
                    Et sauvegarder qui permet de demander à l'utilisateur si il veut enregistrer le résulat.
                    Permet d'obtenir la distance Hamming entre 2 séquences en la parcourant tout le long.


                    Args:
                        - fasta (str): La chaîne de caractères du chemin vers le fichier.
                        - motif (str): La chaîne de caractères qui correspond au motif recherché.
                        - mismatch (int): L'entier qui correspond aux événements de mismatch tolérés.
                        - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

                    Returns:
                        - Aucun.

                    Print:
                        - La position des motifs ressemblant à X mismatch près.

    """
    sequences = lire_fasta.lire(fasta)
    texte = ''
    for sid, seq in sequences.items():
        texte += f"{sid}\n\n"
        for i in range(len(seq)-len(motif)+1):
            if hamming(seq[i:i+len(motif)], motif) <= mismatch:
                texte += f"{i+1}\n{seq[i:i+len(motif)]}\n{motif}\n\n"
        texte += '\n'

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/pattern_approximatif.txt', texte)

    print(texte)


pattern_approximatif('Données/test.txt', 'VHE', 1, True)
