from Module import lire_fasta, code_genetique, enregistrer_donnees


def traduction(fasta: str, sauvegarder: bool = False) -> dict:
    """
            Prend un fichier FASTA.
            Affiche pour chaque séquence du FASTA, la traduction de la séquence.

            Args:
                - fasta (str): La chaîne de caractères qui désigne le chemin du fichier.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                - donnees (dict): Un dictionnaire qui contient {id: peptide}.

            Print:
                - La séquence du peptide.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ""
    code = code_genetique.codon()
    donnees = {}

    for id, sequence in seqs.items():
        prot = ''
        for i in range(0, len(sequence), 3):
            if len(sequence[i:i+3]) == 3 and code[sequence[i:i+3].replace('T', 'U')] == 'STOP':
                break
            elif len(sequence[i:i+3]) == 3:
                prot += code[sequence[i:i+3].replace('T', 'U')]
        texte += f"{id}\n{prot}\n"
        donnees[id] = prot
        print(prot)

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/traduction_résultats_nt_a_ct.txt', texte)

    return donnees


traduction('Résultats/epissage_résultats.txt', True)
