

def table() -> dict:
    """
            Créer un tableau avec le symbole et sa valeur phred.
            Il est basé sur la version ""Illumina 1.8+"".

            Args:
                Aucun.

            Returns:
                - score (dict): Un dictionnaire avec {Symbole: phred}.

    """
    score = {}
    with open('../Module/phred_score_illumina1.8+.txt', 'r') as filin:
        for ligne in filin:
            score[ligne.split()[0]] = int(ligne.split()[1])
    return score
