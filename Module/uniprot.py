import requests


def prots_seqs(id: str) -> str:
    """
        Prend un ID, et renvoie sa séquence sur UNIPROT.

        Args:
            - id (str): La chaîne de caractères correspondant à l'id de la séquence.

        Returns:
            - sequence (str): La chaine de caractères correspondant à la séquence de l'id.

    """
    if '_' in id:
        url = f"https://www.uniprot.org/uniprot/{id[:6]}.fasta"
        response = requests.get(url)
        lines = response.text.split('\n')
        sequence = ''.join(lines[1:])
    else:
        url = f"https://www.uniprot.org/uniprot/{id}.fasta"
        response = requests.get(url)
        lines = response.text.split('\n')
        sequence = ''.join(lines[1:])

    return sequence
