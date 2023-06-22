import math
from In_silico.Module import lire_fastq, table_phred, enregistrer_donnees
import matplotlib.pyplot as plt


def graphique(fastq: str, seuil: int = 2, sauvegarder: bool = False):
    """
            Prend un chemin vers le FASTQ, un seuil qui déterminera la valeur considérée.
            Affiche la moyenne de score pour chaque base, en utilisant l'ensemble des séquences
            présentent dans le FASTQ.

            Args:
                - fastq (str): La chaîne de caractères du chemin vers le fichier.
                - seuil (int): L'entier qui correspond à la limite de valabilité.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                Aucun.

            Print:
                - Le graphique de qualité.

    """
    sequences = lire_fastq.lire(fastq)
    scores = table_phred.table()

    cpt = []
    axe_x = []
    moy_axe_y = []

    for sid, sequence in sequences.items():
        axe_x = max(list(range(1, len(sequence[0])+1)), axe_x, key=len)

        if len(axe_x) > len(moy_axe_y):
            for _ in range(len(axe_x)-len(moy_axe_y)):
                moy_axe_y.append(0)
                cpt.append(0)

        axe_y = []

        for i in range(len(sequence[1])):
            axe_y.append(scores[sequence[1][i]])
            cpt[i] += 1

        for i in range(len(axe_y)):
            moy_axe_y[i] += axe_y[i]

    for i in range(len(moy_axe_y)):
        moy_axe_y[i] /= cpt[i]

    plt.ylim([0, 50])
    plt.xlim([1, len(moy_axe_y)*1.05])
    plt.ylabel('Score de qualité')
    plt.xlabel('Numéro de base')
    plt.title('Graphique représentant le score en fonction de la base')
    plt.plot(axe_x, moy_axe_y, c='black', label='Qualité')
    plt.plot([i for i in range(math.ceil(len(moy_axe_y)*1.05))], [seuil for _ in range(math.ceil(len(moy_axe_y)*1.05))], c='r', label='Seuil')
    plt.gca().legend(('Qualité', 'Seuil'))

    if sauvegarder:
        plt.savefig('Résultats/graphique_qualite_ensemble.pdf')

    plt.show()


def filtre(fastq: str, seuil: int = 20, sauvegarder: bool = True):
    """
                Prend un chemin vers le FASTQ, un seuil qui déterminera la valeur considérée.
                Enregistre dans des fichiers séparés les séquences >= ou < au seuil.

                Args:
                    - fastq (str): La chaîne de caractères du chemin vers le fichier.
                    - seuil (int): L'entier qui correspond à la limite de valabilité.
                    - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

                Returns:
                    Aucun.

                Print:
                    - Les séquences >= au seuil.

    """
    sequences = lire_fastq.lire(fastq)
    scores = table_phred.table()
    texte_inf_seuil = ''
    texte_supp_seuil = ''

    for sid, sequence in sequences.items():
        moyenne = 0

        for lettre in sequence[1]:
            moyenne += scores[lettre]

        moyenne /= len(sequence[1])

        if moyenne < seuil:
            texte_inf_seuil += f'{sid}\n{sequence[0]}\n+\n{sequence[1]}\n'

        else:
            texte_supp_seuil += f'{sid}\n{sequence[0]}\n+\n{sequence[1]}\n'

    if sauvegarder:
        enregistrer_donnees.enregistrer(f"Résultats/sequence_inf_{seuil}_fastq.txt", texte_inf_seuil)
        enregistrer_donnees.enregistrer(f"Résultats/sequence_supp_{seuil}_fastq.txt", texte_supp_seuil)

    print(texte_supp_seuil)


graphique('Données/reads1_fastq.txt', 20, True)
filtre('Données/reads1_fastq.txt', 20, True)
