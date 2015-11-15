import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.db.models import Q

from core.utils import BLAST
from public_interface.models import Sequences


class BLASTNew(BLAST):
    def __init__(self, blast_type, name, sequence, gene_codes):
        """
        :param blast_type: new
        :param name: name given to query sequence
        :param sequence: sequence to blast given by user
        :param gene_codes: list of gene_codes to blast against
        :param mask:
        """
        super(BLASTNew, self).__init__(blast_type, name, sequence, gene_codes)
        self.name = name
        self.sequence = sequence

        if gene_codes:
            gene_codes = [i.gene_code for i in gene_codes]
            gene_codes.sort()
            self.gene_codes = gene_codes
        else:
            self.gene_codes = []

        self.path = os.path.join(self.cwd, 'db',
                                 '_'.join(self.gene_codes) + '_seqs.fas.n*')
        self.db = os.path.join(self.cwd, 'db',
                               '_'.join(self.gene_codes) + '_seqs.fas')

    def save_seqs_to_file(self):
        """Query sequences for each gene from database and save to local disk.

        Sets attribute `self.seq_file` containing necessary sequences from our
        database.

        """
        if self.blast_type == 'new':
            self.seq_file = os.path.join(self.cwd,
                                         'db',
                                         '_'.join(self.gene_codes) + "_seqs.fas",
                                         )
            if self.gene_codes:
                # Taken from http://stackoverflow.com/a/1239602
                Qr = None
                for gene_code in self.gene_codes:
                    q = Q(gene_code=gene_code)
                    if Qr:
                        Qr = Qr | q
                    else:
                        Qr = q
                queryset = Sequences.objects.filter(Qr)
            else:
                queryset = Sequences.objects.all()

            my_records = []
            for i in queryset:
                item_id = i.code_id + '|' + i.gene_code
                seq = self.strip_question_marks(i.sequences)
                if seq != '':
                    seq_record = SeqRecord(Seq(seq), id=item_id)
                    my_records.append(seq_record)
            SeqIO.write(my_records, self.seq_file, "fasta")

    def save_query_to_file(self):
        this_id = self.name
        seq = self.strip_question_marks(self.sequence)

        if seq != '':
            seq_record = SeqRecord(Seq(seq),
                                   id=this_id)
            SeqIO.write(seq_record, self.query_file, "fasta")
