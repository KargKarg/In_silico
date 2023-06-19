def tableau_masse() -> dict:
    """
        Créer un tableau avec l'AA et sa monoisotopic masse.

        Args:
            Aucun.

        Returns:
            - tableau (dict): Un dictionnaire avec {AA: monoisotopic masse}.

    """
    tableau = {}
    with open('../Module/table_poids_aa.txt', 'r') as filin:
        for ligne in filin:
            tableau[ligne.split()[0]] = float(ligne.split()[1])
    return tableau
