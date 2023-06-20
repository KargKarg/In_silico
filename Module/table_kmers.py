import itertools


def kmers(k: int) -> dict:
    """
            Prend un entier et renvoie un tableau avec K-mers et occurence.
            Les occurences sont initialisés à 0.

            Args:
                - k (int): L'entier définissant la taille des k-mers.

            Returns:
                - tableau (dict): Un dictionnaire avec {kmers: 0}.

    """
    bases = ['A', 'C', 'G', 'T']
    tableau = {''.join(combi): 0 for combi in itertools.product(bases, repeat=k)}
    return tableau
