from In_silico.Module import lire_fasta, enregistrer_donnees
import brin_comp


def distance_hamming(sequence1: str, sequence2: str) -> tuple:
    """
                Renvoie la distance Hamming des deux séquences, la distance est incrémentée lorsque
                sequence1[i] != sequence2[i].

                Args:
                    - sequence1 (str): La chaîne de caractères de la première séquence.
                    - sequence2 (str): La chaîne de caractères de la deuxième séquence.

                Returns:
                    - La distance Hamming et l'indice de la dernière mutation observée.

                Print:
                    - Aucun.

    """
    distance = 0
    indice = 0
    for i in range(min(len(sequence1), len(sequence2))):
        if sequence1[i] != sequence2[i]:
            distance += 1
            indice = i
    return distance, indice


def correction(fasta: str, sauvegarder: bool = False):
    """
                Affiche les reads corrigés, lorsque la read ne possède qu'une erreur et qu'elle n'est présente
                qu'une seule fois dans le fichier.

                Args:
                    - fasta (str): La chaîne de caractères du chemin du fichier FASTA.
                    - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

                Returns:
                    - Aucun.

                Print:
                    - Affiche la correction avec read_erronnée->read_corrigée.
    """
    sequences = lire_fasta.lire(fasta)
    sequences_compl = brin_comp.brin_complementaire(fasta).values()
    sequences_erronees = []
    sequences_corrigees = []
    texte = ''

    for sequence_indirect in sequences_compl:
        sequence_direct = sequence_indirect[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

        if sequence_indirect not in sequences.values() and list(sequences.values()).count(sequence_direct) < 2:
            sequences_erronees.append(sequence_direct)

    for erreur in sequences_erronees:
        for reference in sequences.values():

            reference_indirect = reference[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

            if erreur != reference and distance_hamming(erreur, reference)[0] == 1 and reference not in sequences_erronees:
                _, indice = distance_hamming(erreur, reference)
                sequences_corrigees.append(f"{erreur[:indice]}{reference[indice]}{erreur[indice+1:]}")
                break

            elif erreur != reference and distance_hamming(erreur, reference_indirect)[0] == 1 and reference not in sequences_erronees:
                _, indice = distance_hamming(erreur, reference_indirect)
                sequences_corrigees.append(f"{erreur[:indice]}{reference_indirect[indice]}{erreur[indice+1:]}")
                break

    for i in range(len(sequences_erronees)):
        texte += f"{sequences_erronees[i]}->{sequences_corrigees[i]}\n"

    if sauvegarder:
        enregistrer_donnees.enregistrer('Résultats/correction_reads.txt', texte)

    print(texte)


correction('Données/test.txt', True)
