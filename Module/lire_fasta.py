
def lire(fasta: str) -> dict:
    """
        Prend un chemin vers le fichier FASTA.
        Créer un tableau avec id et séquence.

        Args:
            - fasta (str): La chaîne de caractères du chemin vers le fichier.

        Returns:
            - data (dict): Un dictionnaire avec {id de la séquence: séquence}.

    """
    data = {}

    with open(fasta, 'r') as filin:
        id = ''
        seq = ''
        for ligne in filin:
            if ligne[0] == '>':
                data[id] = seq
                try:
                    id = ligne.split('|')[1]
                except IndexError:
                    id = ligne.split()[0]
                seq = ''
            else:
                seq += ligne.replace('\n', '')

    data[id] = seq
    del data['']
    return data
