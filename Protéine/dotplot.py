from In_silico.Module import lire_fasta, Levenshtein
import matplotlib.pyplot as plt


def dotplot(fasta: str, fenetre: int, pourcentage_identite: float = 100, sauvegarder: bool = False):
    """
            Prend un chemin vers le fichier FASTA.
            Une fenêtre qui glissera le long de la séquence.
            Un pourcentage d'identité, si la valeur de la distance Levenshtein >=, on considére une identité.
            Et sauvegarder qui permet de demander à l'utilisateur si il veut enregistrer le résulat.


            Args:
                - fasta (str): La chaîne de caractères du chemin vers le fichier.
                - fenetre (int): L'entier qui donnera la taille de la fenêtre.
                - pourcentage_identite (float): Donne la valeur seuil à dépasser.
                - sauvegarder (bool): Si True, le fichier est sauvegardé dans le répertoire "Résultats".

            Returns:
                - Aucun.

            Print:
                - Le DotPlot entre les deux séquences.

    """
    nomseq1, nomseq2 = lire_fasta.lire(fasta).keys()
    seq1, seq2 = lire_fasta.lire(fasta).values()

    axe_x = []
    axe_y = []

    for i in range(len(seq1)-fenetre+1):
        for j in range(len(seq2)-fenetre+1):
            if Levenshtein.pourcentage_ressemblance(seq1[i:i+fenetre], seq2[j:j+fenetre]) >= pourcentage_identite:
                print(seq1[i:i+fenetre], seq2[j:j+fenetre])
                axe_x.append(i+1)
                axe_y.append(j+1)

    plt.scatter(axe_x, axe_y, marker='o', s=5, c='m')
    plt.xlim(0, len(seq1))
    plt.ylim(0, len(seq2))
    plt.title(f"Dotplot entre {nomseq1[1:]} et {nomseq2[1:]}\n"
              f"avec fenetre = {fenetre} et pourcentage pour identite = {pourcentage_identite}%")
    plt.xlabel(nomseq1[1:])
    plt.ylabel(nomseq2[1:])

    if sauvegarder:
        plt.savefig('Résultats/dotplot.pdf')

    plt.show()


dotplot("Données/test.txt", 4, 75, True)
