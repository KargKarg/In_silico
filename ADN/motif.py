from Module import lire_fasta, enregistrer_donnees
import re


def recherche_motif(fasta: str, motif: str, sauvegarder: bool = False):
    """
        Prend un chemin vers le fichier FASTA, un motif d'intérêt.
        Affiche pour chaque séquence de ce FASTA, la position du motif.

        A noter qu'un REGEX peut être utilisé.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - motif (str): La chaîne de caractère qui correspond au motif recherché.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Le résultat de la recherche du motif.

    """
    seqs = lire_fasta.lire(fasta)
    texte = ""
    regex = re.compile(motif)
    for id, sequence in seqs.items():
        texte += f'{id}\n'
        occurences = re.finditer(regex, sequence)
        for occurence in occurences:
            texte += f'{occurence.start()+1} '
        texte += '\n'
    if sauvegarder:
        enregistrer_donnees.enregistrer(f"Résultats/motif_résultats_{motif}.txt", texte)
    print(texte)


recherche_motif('Résultats/brin_complémentaire_résultats_5a3.txt', 'TTGACA')
