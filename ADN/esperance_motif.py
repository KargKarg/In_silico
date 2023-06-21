from In_silico.Module import lire_fasta, enregistrer_donnees
import math


def proportion_gc(sequence: str) -> float:
    """
            Prend une séquence et donne sa proportion en GC.

            Args:
                - sequence (str): La chaîne de caractères correspondant à la séquence.

            Returns:
                - (float): Proportion en GC.

            Print:
                - Aucun.

    """
    return (sequence.count('G') + sequence.count('C'))/len(sequence)


def probabilite(motif: str, contenu_gc: float) -> float:
    """
            Prend un motif et calcul la probabilité que ce motif apparaisse sur une séquence qui possède ce contenu_gc.

            Args:
                - motif (str): La chaîne de caractères correspondant au motif.
                - contenu_gc (bool): La proportion en GC d'une séquence.

            Returns:
                - (float): La probabilité d'obtenir le motif.

            Print:
                - Aucun.

    """
    return 1 * math.pow((1-contenu_gc)/2, motif.count('A')
                        + motif.count('T')) * math.pow(contenu_gc/2, motif.count('C')+motif.count('G'))


def esperance(fasta: str, motif: str, sauvegarder: bool = False):
    """
            Prend un chemin vers le fichier FASTA et un motif pour nous donner l'espérance du motif dans la séquence.

            Args:
                - fasta (str): La chaîne de caractères du chemin vers le fichier.
                - motif (str): La chaîne de caractères correspondant au motif.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                - Aucun.

            Print:
                - L'espérance pour chaque séquence dans le FASTA avec le motif affiché.

    """
    sequences = lire_fasta.lire(fasta)
    texte = ''
    for sid, sequence in sequences.items():
        content = proportion_gc(sequence)
        proba_motif = probabilite(motif, content)
        esper = proba_motif*(len(sequence)-len(motif)+1)
        texte += f"{sid} {round(esper, 3)} {motif}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/esperance_motif.txt', texte)

    print(texte)


esperance('Données/test.txt', 'GGCC', True)
