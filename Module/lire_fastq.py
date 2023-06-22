
def lire(fastq: str) -> dict:
    """
        Prend un chemin vers le fichier FASTQ.
        Créer un tableau avec id et [séquence, séquence qualité].

        Args:
            - fastq (str): La chaîne de caractères du chemin vers le fichier.

        Returns:
            - sequences (dict): Un dictionnaire avec {id de la séquence: [séquence, séquence qualité]}.

    """
    sequences = {}

    with open(fastq, 'r') as filin:
        for ligne in filin:
            if ligne.replace('\n', '') not in sequences.keys():
                sequences[ligne.replace('\n', '')] = ['', '']
            data = []
            for _ in range(3):
                data.append(filin.readline().replace('\n', ''))
            sequences[ligne.replace('\n', '')][0] += data[0]
            sequences[ligne.replace('\n', '')][1] += data[2]

    return sequences
