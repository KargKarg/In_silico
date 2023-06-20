from In_silico.Module import lire_fasta, enregistrer_donnees


def transitions_transversions(fasta: str, sauvegarder: bool = False):
    """
            Prend un chemin vers le fichier FASTA.
            Affiche les positions des transitions/transversions et donne ce même ratio.

            Args:
                - fasta (str): La chaîne de caractères du chemin vers le fichier.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                Aucun.

            Print:
                - Le résultat des événements transitions/transversions.

    """
    sequence1, sequence2 = lire_fasta.lire(fasta).values()
    pattern = [('A', 'G'), ('C', 'T')]
    transitions = 0
    transversions = 0
    texte = ''

    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            if (sequence1[i], sequence2[i]) in pattern or (sequence2[i], sequence1[i]) in pattern:
                transitions += 1
                texte += f"transition {i+1}pb\n"
            else:
                transversions += 1
                texte += f"transversion {i+1}pb\n"

    texte += f"Ratio de transitions/transversions : {transitions/transversions}"

    if sauvegarder:
        enregistrer_donnees.enregistrer("Résultats/transitions_transversions.txt", texte)

    print(texte)


transitions_transversions('Données/test.txt', True)
