import glob
import os
import time

from django.core.management import call_command
from django.test import TestCase
from django.conf import settings

from .utils import BLAST
from public_interface.models import Vouchers
from public_interface.models import Sequences


class BlastLocalTest(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        blast_type = 'local'
        voucher_code = 'CP100-10'
        gene_code = 'COI'
        self.blast = BLAST(blast_type, voucher_code, gene_code, test=True)
        self.seq_file = ''

    def remove_blast_data_files(self):
        path = os.path.join(settings.BASE_DIR,
                            '..',
                            'blast_local',
                            'db',
                            '*',
                            )
        for file in glob.glob(path):
            if not file.endswith('py') and os.path.isfile(file):
                os.remove(file)

    def test_have_blast_db(self):
        self.remove_blast_data_files()
        result = self.blast.have_blast_db()
        self.assertEqual(False, result)

    def test_have_blast_db_true(self):
        self.blast.save_seqs_to_file()
        self.blast.create_blast_db()

        result = self.blast.have_blast_db()
        self.assertEqual(True, result)
        self.remove_blast_data_files()

    def test_save_seqs_to_file(self):
        self.blast.save_seqs_to_file()
        self.seq_file = os.path.join(settings.BASE_DIR,
                                     '..',
                                     'blast_local',
                                     'db',
                                     'COI_seqs.fas',
                                     )
        result = os.path.isfile(self.seq_file)
        self.assertTrue(result)

    def test_create_blast_db(self):
        blast_type = 'local'
        voucher_code = 'CP100-10'
        gene_code = 'COI'
        blast = BLAST(blast_type, voucher_code, gene_code, test=True, mask=True)
        blast.save_seqs_to_file()
        blast.create_blast_db()
        files = glob.glob(
            os.path.join(settings.BASE_DIR,
                         '..',
                         'blast_local',
                         'db',
                         'COI_seqs.fas.n*')
        )
        self.assertTrue(len(files) > 0)
        self.remove_blast_data_files()

    def test_create_blast_db_unmasked(self):
        blast_type = 'local'
        voucher_code = 'CP100-10'
        gene_code = 'COI'
        blast = BLAST(blast_type, voucher_code, gene_code, test=True, mask=False)
        blast.save_seqs_to_file()
        blast.create_blast_db()
        files = glob.glob(
            os.path.join(settings.BASE_DIR,
                         '..',
                         'blast_local',
                         'db',
                         'COI_seqs.fas.n*')
        )
        self.assertTrue(len(files) > 0)
        self.remove_blast_data_files()

    def test_is_blast_db_up_to_date(self):
        self.blast.save_seqs_to_file()
        self.blast.create_blast_db()
        result = self.blast.is_blast_db_up_to_date()
        self.assertTrue(result)

    """
    def test_is_blast_db_up_to_date_false(self):
        self.blast.save_seqs_to_file()
        self.blast.create_blast_db()

        time.sleep(5)
        b = Vouchers.objects.get(code='CP100-10')
        c = Sequences.objects.get(code=b, gene_code='COI')
        print(">>>>c.time_edited", c.time_edited)

        Sequences.objects.filter(code=b, gene_code='COI').update(
            notes="a random update for this record."
        )
        c = Sequences.objects.get(code=b, gene_code='COI')
        print(">>>>c.time_edited", c.time_edited)
        result = self.blast.is_blast_db_up_to_date()
        self.assertFalse(result)
    """
