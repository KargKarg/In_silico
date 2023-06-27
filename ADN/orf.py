from In_silico.Module import lire_fasta, code_genetique, enregistrer_donnees
import brin_comp


def chercher_orf(fasta: str, sauvergarder: bool = False) -> None:
    """
        Prend un chemin vers le fichier FASTA.
        Affiche tous les ORFs des séquences dans le FASTA:
            - brin direct ou indirect.
            - position [début, fin].
            - cadre de lecture.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.
            - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

        Returns:
            Aucun.

        Print:
            - Le résultat de la recherche des ORFs.

    """
    seqs = lire_fasta.lire(fasta)
    seqs_comp = brin_comp.brin_complementaire(fasta)
    code = code_genetique.codon()
    texte = ""
    orf = 0
    for (id, sequence), (_, sequence_comp) in zip(seqs.items(), seqs_comp.items()):
        for i in range(len(sequence)):
            if sequence[i:i+3] == 'ATG':
                orf += 1
                for k in range(i+3, len(sequence), 3):
                    if len(sequence[k:k+3]) == 3 and code[sequence[k:k+3].replace('T', 'U')] == 'STOP':
                        texte += f"{id}_ORF{orf} brin direct [{i+1};{k+3}] cadre : {i%3}\n{sequence[i:k+3]}\n"
                        break
            if sequence_comp[i:i+3] == 'ATG':
                orf += 1
                for k in range(i+3, len(sequence_comp), 3):
                    if len(sequence_comp[k:k+3]) == 3 and code[sequence_comp[k:k+3].replace('T', 'U')] == 'STOP':
                        texte += f"{id}_ORF{orf} brin indirect " \
                                 f"[{len(sequence)-i};{len(sequence)-k-2}] cadre : {i%3}\n{sequence_comp[i:k+3]}\n"
                        break
        texte += '\n'

    if sauvergarder:
        enregistrer_donnees.enregistrer("Résultats/orf_résultats.txt", texte)

    print(texte)


chercher_orf('Données/seq.txt', True)
