
def codon() -> dict:
    """
        Permet d'avoir un tableau avec pour chaque codon, son AA associé.

        Args:
            Aucun.

        Returns:
            - code_g (dict):  Un dictionnaire avec {Codon: AA}.

    """
    code_g = {}
    with open('../Module/code_genetique.txt', 'r') as filin:
        for ligne in filin:
            ligne = ligne.split()
            code_g[ligne[0]] = ligne[1]
    return code_g
