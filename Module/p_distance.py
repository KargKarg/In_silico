

def calcul(sequence1: str, sequence2: str) -> float:
    """
                    Prend deux séquences, et calcul la ""p distance"" entre les deux séquences.
                    La ""la pidstance"" correspond au nombre d'élément différent divisé par la taille de la séquence.

                    Args:
                        - chaine1 (string): Première chaîne.
                        - chaine2 (string): Deuxième chaîne.

                    Returns:
                        - (float): p distance.

                    Print:
                        - Aucun.

    """
    difference = 0

    for i in range(min(len(sequence1), len(sequence2))):
        if sequence1[i] != sequence2[i]:
            difference += 1

    return round(difference/min(len(sequence1), len(sequence2)), 5)
