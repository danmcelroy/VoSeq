import uuid
import os


class DatasetHandler(object):
    def __init__(self, dataset_str, file_format):
        self.dataset_str = dataset_str
        self.file_format = file_format
        self.warnings = []
        self.errors = []

        self.cwd = os.path.dirname(__file__)
        self.guid = self.make_guid()
        self.dataset_file = os.path.join(self.cwd,
                                         'dataset_files',
                                         self.file_format + '_' + self.guid + '.txt',
                                         )
        self.aa_dataset_file = os.path.join(self.cwd,
                                            'dataset_files',
                                            self.file_format + '_aa_' + self.guid + '.txt',
                                            )

    def save_dataset_to_file(self):
        with open(self.dataset_file, 'w') as handle:
            handle.write(self.dataset_str)

    def make_guid(self):
        return uuid.uuid4().hex
