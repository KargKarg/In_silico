from Module import lire_fasta, enregistrer_donnees
import brin_comp

def semi_palindrome(fasta: str, sauvegarder: bool = True):
    """
        Prend un chemin vers le fichier FASTA.
        Affiche pour chaque séquence du FASTA, les motifs semi-palindromique et leurs positions.

        Peut servir à détécter des sites de restriction.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Les positions et motifs des sequences semi-palindromiques.

    """
    seqs = lire_fasta.lire(fasta)
    seqs_comp = brin_comp.brin_complementaire(fasta)
    texte = ""
    for id, sequence in seqs.items():
        texte += f'{id}\n'
        for i in range(len(sequence)):
            for j in range(i+4, i+14, 2):
                if sequence[i+(j-i)//2:j] == seqs_comp[id][::-1][i:i+(j-i)//2][::-1]:
                    texte += f"{i + 1} {sequence[i:i + (j - i) // 2]} {sequence[i + (j - i) // 2:j]} {j} \n"
    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/semi_palindrome_résultats.txt', texte)

    print(texte)


semi_palindrome('Données/seq.txt')
