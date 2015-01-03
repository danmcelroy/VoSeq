import os
import re
import subprocess

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from public_interface.models import Sequences


class BLAST(object):
    """
    Class to handle duties related to local blast against sequences of one gene,
    and full blast against all sequences in our database.
    """
    def __init__(self, blast_type, voucher_code, gene_code):
        """
        Type of blast to do: local, full, remote

        :param blast_type: local, full, remote.
        :param voucher_code:
        :param gene_code:
        """
        self.blast_type = blast_type
        self.voucher_code = voucher_code
        self.gene_code = gene_code
        self.mask = True
        self.cwd = os.path.dirname(__file__)

    def have_blast_db(self):
        """
        Finds out whether we already have a blast db with our sequences.

        :return: True or False
        """
        pass

    def is_blast_db_up_to_date(self):
        """
        Finds out whether our blast db contains all our sequences. In other
        words, it finds out whether there are sequences in our postgres db with
        time_created or time_edited more recent than our blast db files.

        :return:
        """
        pass

    def save_seqs_to_file(self):
        """
        Query sequences for each gene from our database and save them to local
        disk.

        Sets attribute `self.seq_file` containing necessary sequences from our
        database.
        """
        if self.blast_type == 'local':
            self.seq_file = os.path.join(self.cwd,
                                         'db',
                                         self.gene_code + "_seqs.fas",
                                         )
            queryset = Sequences.objects.all().filter(gene_code=self.gene_code)

            my_records = []
            for i in queryset:
                id = i.code_id + '|' + i.gene_code
                seq = self.strip_question_marks(i.sequences)
                seq_record = SeqRecord(Seq(seq),
                                       id=id)
                my_records.append(seq_record)
            SeqIO.write(my_records, self.seq_file, "fasta")

    def strip_question_marks(self, seq):
        seq = re.sub('^\?+', '', seq)
        seq = re.sub('\?+$', '', seq)
        return seq

    def create_blast_db(self):
        """
        Creates a BLAST database from our sequences file in FASTA format.
        Optionally eliminates low-complexity regions from the sequences.

        :return:
        """
        if self.mask is True:
            command = 'dustmasker -in ' + self.seq_file + ' -infmt fasta '
            command += '-outfmt maskinfo_asn1_bin -out ' + self.seq_file + '_dust.asnb'
            subprocess.check_output(command, shell=True)  # identifying low-complexity regions.

            command = 'makeblastdb -in ' + self.seq_file + ' -input_type fasta -dbtype nucl '
            command += '-mask_data ' + self.seq_file + '_dust.asnb '
            command += '-out ' + self.seq_file + ' -title "Whole Genome without low-complexity regions"'
            print("creating database...")
            subprocess.check_output(command, shell=True)  # Overwriting the genome file.
        else:
            command = 'makeblastdb -in ' + self.seq_file + ' -input_type fasta -dbtype nucl '
            command += '-out ' + self.seq_file + ' -title "Whole Genome unmasked"'
            print("creating database...")
            subprocess.check_output(command, shell=True)

    def do_blast(self):
        blastn_cline = NcbiblastnCommandline(query=self.query, db=self.db,
                                             evalue=0.001, outfmt=5, out="opuntia.xml")
        blastn_cline()
