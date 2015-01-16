import os
import uuid

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from blast_local.utils import BLAST
from public_interface.models import Sequences


class BLASTNew(BLAST):
    def __init__(self, blast_type, name, sequence, gene_codes, mask=None):
        """
        :param blast_type: new
        :param name: name given to query sequence
        :param sequence: sequence to blast given by user
        :param gene_codes: list of gene_codes to blast against
        :param mask:
        """
        self.e_value = 0.001
        self.blast_type = blast_type
        self.name = name
        self.sequence = sequence
        self.gene_codes = gene_codes
        self.cwd = os.path.dirname(__file__)
        self.seq_file = ""

        if mask is not False:
            self.mask = True
        else:
            self.mask = False

        self.path = os.path.join(self.cwd,
                                 'db',
                                 '_'.join(self.gene_codes) + '_seqs.fas.n*',
                                 )
        self.db = os.path.join(self.cwd,
                               'db',
                               '_'.join(self.gene_codes) + '_seqs.fas',
                               )
        self.query_file = os.path.join(self.cwd,
                                       'db',
                                       'query_' + uuid.uuid4().hex + '.fas',
                                       )
        self.output_file = os.path.join(self.cwd,
                                        'db',
                                        'output_' + uuid.uuid4().hex + '.xml',
                                        )
