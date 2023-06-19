from Module import lire_fasta, enregistrer_donnees


def motifs_partages(fasta: str, sauvergarder: bool = False):
    """
        Prend un chemin vers le fichier FASTA.
        Affiche le plus long motif qui est présent dans toutes les séquences du FASTA.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Le motif partagé le plus long.

    """
    seqs = lire_fasta.lire(fasta)
    motif = ''
    taille = len(motif)
    fid, fsequence = iter(seqs.items()).__next__()
    for i in range(len(fsequence)):
        contenu = True
        for j in range(len(fsequence)+1, i+1, -1):
            for sid, ssequence in seqs.items():
                if fid != sid and ssequence.find(fsequence[i:j]) == -1 and contenu:
                    contenu = False
                    break
            if contenu and len(fsequence[i:j]) > taille:
                motif = fsequence[i:j]
                taille = len(fsequence[i:j])
                break
        if contenu and len(fsequence)-i < taille:
            break
    if sauvergarder:
        enregistrer_donnees.enregistrer('Résultats/motif_partage_résultats.txt', motif)

    print(motif)


motifs_partages('Données/seq.txt', True)
