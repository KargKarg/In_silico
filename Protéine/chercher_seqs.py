from Module import uniprot, enregistrer_donnees


def chercher_seqs(all_ids: list, sauvegarder: bool = False) -> dict:
    """
        Prend une liste d'IDs, et renvoie un tableau avec id et séquence associée.

        Args:
            - all_ids (list): Liste contenant l'ensemble des IDs.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            - seqs (dict): Un dictionnaire contenant {id: sequence}

    """
    texte = ''
    seqs = {}
    for id in all_ids:
        seq = uniprot.prots_seqs(id)
        seqs[id] = seq
        texte += f">{id}\n{seq}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/sequences_résultats_uniprot.txt', texte)

    return seqs


chercher_seqs(['A2Z669', 'B5ZC00', 'P07204_TRBM_HUMAN', 'P20840_SAG1_YEAST'], True)
