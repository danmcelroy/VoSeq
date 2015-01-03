import os

from django.test import TestCase
from django.conf import settings

from .utils import BLAST


class BlastLocalTest(TestCase):
    def setUp(self):
        blast_type = 'local'
        voucher_code = 'CP100-10'
        gene_code = 'COI'
        self.blast = BLAST(blast_type, voucher_code, gene_code)

    def test_have_blast_db(self):
        result = self.blast.have_blast_db()
        self.assertEqual(False, result)

    def test_have_blast_db_true(self):
        self.blast.create_blast_db()

        result = self.blast.have_blast_db()
        self.assertEqual(True, result)

    def test_save_seqs_to_file(self):
        self.blast.save_seqs_to_file()
        expected_file = os.path.join(settings.BASE_DIR,
                                     '..',
                                     'blast_local',
                                     'db',
                                     'COI_seqs.fas',
                                     )
        result = os.path.isfile(expected_file)
        self.assertTrue(result)
