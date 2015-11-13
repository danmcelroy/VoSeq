from io import StringIO
import os
import unittest
from unittest.mock import patch

from django.conf import settings
from django.test import TestCase
from django.core.management import call_command

from blast_ncbi.utils import BLASTNcbi


TEST_PATH = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(TEST_PATH, "CP100-10_COI-begin.xml"), "r") as handle:
    ncbi_return_handle = StringIO(handle.read())


class TestNcbiBlast(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump2.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        voucher_code = 'CP100-10'
        gene_code = 'COI-begin'
        self.blast = BLASTNcbi(voucher_code, gene_code)

    @patch("Bio.Blast.NCBIWWW.qblast", return_value=ncbi_return_handle)
    def test_blast_with_accession_number_in_header(self, mock_qblast):
        self.blast.save_query_to_file()
        self.blast.do_blast()
        result = self.blast.parse_blast_output()
        self.blast.delete_query_output_files()
        self.assertTrue(len(result) > 0)

    @unittest.skipIf(settings.TRAVIS is True,
                     'Testing using BLASTNcbi fails due to network problems')
    def test_index(self):
        response = self.client.get('/blast_ncbi/CP100-10/COI/')
        self.assertEqual(200, response.status_code)
