class CreateDataset(object):
    """
    Accept form input to create a dataset in several formats, codon positions,
    for list of codes and genes.

    Attributes:
        ``dataset_str``: output dataset to pass to users.

    """
    def __init__(self):
        self.dataset_str = self.create_dataset()

    def create_dataset(self):
        return "hola"
