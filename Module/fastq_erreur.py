import math
from In_silico.Module import table_phred


def proba(symbole: str) -> float:
    """
                Calcul la probabilité d'erreur avec 10^(-P/10).

                Args:
                    - symbole (str): Correspond au symbole associé à une base.

                Returns:
                    - (float): La probabilité d'erreur.

    """
    tableau = table_phred.table()
    return math.pow(10, (-tableau[symbole]/10))
