
def enregistrer(path: str, texte: str) -> None:
    """
        Prend un chemin vers un fichier.
        Ecrit dans le fichier le texte passé.

        Args:
            - path (str): La chaîne de caractères du chemin vers le fichier.
            - texte (str): La chaîne de caractères devant figurer dans le fichier.

        Returns:
            Aucun.

    """
    with open(path, 'w') as filout:
        filout.write(texte)
