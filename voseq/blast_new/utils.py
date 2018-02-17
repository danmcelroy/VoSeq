import logging
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.db.models import Q, QuerySet

from core.utils import BLAST
from public_interface.models import Sequences


log = logging.getLogger(__name__)


class BLASTNew(BLAST):
    def __init__(self, blast_type: str, name: str, sequence: str,
                 gene_codes: QuerySet) -> None:
        """
        :param blast_type: new
        :param name: name given to query sequence
        :param sequence: sequence to blast given by user
        :param gene_codes: queryset of Genes to blast against
        """
        super(BLASTNew, self).__init__(blast_type, name, sequence, gene_codes)
        self.name = name
        self.sequence = sequence

        log.debug("Will do blastnew of %s", str(gene_codes))
        if gene_codes:
            self.genes = gene_codes.order_by(
                "gene_code").values_list("gene_code", flat=True)
        else:
            self.genes = []

        self.path = os.path.join(self.cwd, 'db', '_'.join(self.genes) + '_seqs.fas.n*')
        self.db = os.path.join(self.cwd, 'db', '_'.join(self.genes) + '_seqs.fas')

    def save_seqs_to_file(self):
        """Query sequences for each gene from database and save to local disk.

        Sets attribute `self.seq_file` containing necessary sequences from our
        database.

        """
        if self.blast_type == 'new':
            self.seq_file = os.path.join(
                self.cwd,
                'db',
                '_'.join(self.genes) + "_seqs.fas",
            )
            if self.genes:
                # Taken from http://stackoverflow.com/a/1239602
                Qr = None
                for gene_code in self.genes:
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
            seq_record = SeqRecord(Seq(seq), id=this_id)
            SeqIO.write(seq_record, self.query_file, "fasta")
