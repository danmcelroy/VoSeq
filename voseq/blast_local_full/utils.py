import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import BLAST
from public_interface.models import Sequences


class BLASTFull(BLAST):
    """Handles local blast against all sequences in our database.

    The database is `masked` by default, to eliminate low-complexity regions
    from the sequences.

    Use `mask=False` to create unmasked blast databases.

    """
    def __init__(self, *args, **kwargs):
        super(BLASTFull, self).__init__(*args, **kwargs)
        self.path = os.path.join(self.cwd, 'db', 'full_db_seqs.fas.n*')
        self.db = os.path.join(self.cwd, 'db', 'full_db_seqs.fas')

    def save_seqs_to_file(self):
        """Query all sequences from our database and save them to local disk.

        Sets attribute `self.seq_file` containing necessary sequences from our
        database.

        """
        if self.blast_type == 'full':
            self.seq_file = os.path.join(self.cwd, 'db', 'full_db_seqs.fas')
            queryset = Sequences.objects.all()

            my_records = []
            for i in queryset:
                item_id = i.code_id + '|' + i.gene_code
                seq = self.strip_question_marks(i.sequences)
                if seq != '':
                    seq_record = SeqRecord(Seq(seq), id=item_id)
                    my_records.append(seq_record)
            SeqIO.write(my_records, self.seq_file, "fasta")
