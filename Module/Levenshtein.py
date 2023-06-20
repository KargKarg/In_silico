import editdistance


def pourcentage_ressemblance(chaine1: str, chaine2: str) -> float:
    """
            Permet d'obtenir le pourcentage de ressemblance des deux chaînes de caractères,
            grâce à la méthode de Levenshtein.

            Args:
                - chaine1 (string): Première chaîne.
                - chaine2 (string): Deuxième chaîne.

            Returns:
                - pourcentage (flaot):  Le pourcentage de ressemblance des deux chaînes.

    """
    distance = editdistance.eval(chaine1, chaine2)
    taille_max = max(len(chaine1), len(chaine2))
    pourcentage = (1 - distance / taille_max) * 100
    return pourcentage

