import lire_fasta
import enregistrer_donnees


def concatener(fastas: list):
    """
            Prend une liste de chemins vers les différents fichiers FASTA.
            Et crée un unique FASTA avec toutes les séquences.

            Args:
                - fastas (list): La liste contenant les paths des FASTA.

            Returns:
                - Aucun.

    """
    texte = ''
    for path in fastas:
        sequences = lire_fasta.lire(path)
        for sid, sequence in sequences.items():
            texte += f"{sid}\n{sequence}\n"
    enregistrer_donnees.enregistrer('concatener.txt', texte)


concatener(['gene1.txt', 'gene2.txt', 'gene3.txt'])
