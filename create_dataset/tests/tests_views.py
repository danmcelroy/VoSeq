import re

from django.conf import settings
from django.core.management import call_command
from django.db import connection
from django.test import TestCase
from django.test.client import Client

from public_interface.models import Genes, Sequences, Vouchers
from django.contrib.auth.models import User


class CreateDatasetViewsTest(TestCase):
    def setUp(self):
        with connection.cursor() as cursor:
            cursor.execute("alter sequence public_interface_genesets_id_seq restart with 1")
            cursor.execute("alter sequence public_interface_taxonsets_id_seq restart with 1")
            cursor.execute("alter sequence public_interface_sequences_id_seq restart with 1")
        args = []
        opts = {'dumpfile': settings.MEDIA_ROOT + 'test_data.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        g1 = Genes.objects.get(gene_code='COI-begin')
        g2 = Genes.objects.get(gene_code='ef1a')
        self.cleaned_data = {
            'gene_codes': [g1, g2],
            'taxonset': None,
            'voucher_codes': 'CP100-10\r\nCP100-11',
            'geneset': None,
            'outgroup': '',
        }

        self.c = Client()
        self.user = User.objects.get(username='admin')
        self.user.set_password('pass')
        self.user.save()
        self.maxDiff = None

    def test_view_index(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        res = self.c.get('/create_dataset/')
        self.assertEqual(200, res.status_code)

    def test_view_result(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        res = self.c.post('/create_dataset/results/',
                          {
                              'voucher_codes': 'CP100-10',
                              'gene_codes': [],
                              'geneset': 1,
                              'taxonset': 1,
                              'introns': 'YES',
                              'positions': 'ALL',
                              'translations': False,
                              'partition_by_positions': 'ONE',
                              'file_format': 'FASTA',
                              'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                              'degen_translations': 'NORMAL',
                              'exclude': 'YES',
                              'aminoacids': False,
                              'special': False,
                              'outgroup': '',
                          }
                          )
        self.assertEqual(200, res.status_code)
        self.assertEqual('', res)

    def test_view_result_invalid_form(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        res = self.c.post('/create_dataset/results/',
                          {
                              'voucher_codes': '',
                              'gene_codes': [],
                              'geneset': 1,
                              'taxonset': 1,
                          }
                          )
        self.assertEqual(200, res.status_code)

    def test_view_result_get(self):
        res = self.c.get('/create_dataset/results/')
        self.assertEqual(302, res.status_code)

    def test_view_getting_file(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        res = self.c.post('/create_dataset/results/',
                          {
                              'voucher_codes': 'CP100-10',
                              'gene_codes': [],
                              'geneset': 1,
                              'taxonset': 1,
                              'introns': 'YES',
                              'positions': 'ALL',
                              'translations': False,
                              'partition_by_positions': 'by gene',
                              'file_format': 'FASTA',
                              'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                              'degen_translations': 'normal',
                              'exclude': 'YES',
                              'aminoacids': False,
                              'special': False,
                              'outgroup': '',
                          }
                          )
        html_page = res.content.decode('utf-8')
        file_name = re.search('FASTA_\w+\.txt', html_page).group()
        file_content = self.c.get('/create_dataset/results/' + file_name, follow=True)
        expected = ">CP100_10_Aus_aus\nACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACG"
        self.assertTrue(expected in file_content.content.decode('utf-8'))

    def test_view_attemp_to_create_dataset_aa_with_bad_codon(self):
        """Test when trying to translate 'N--' codon. Should translate to X"""
        v = Vouchers.objects.get(code="CP100-10")
        seq = Sequences.objects.get(code=v, gene_code="COI-begin")
        seq.sequences = "TCAN--CGTCCC"
        seq.save()

        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        res = self.c.post('/create_dataset/results/',
                          {
                              'voucher_codes': 'CP100-10',
                              'gene_codes': [],
                              'geneset': 1,
                              'taxonset': 1,
                              'introns': 'YES',
                              'positions': 'ALL',
                              'translations': False,
                              'partition_by_positions': 'by gene',
                              'file_format': 'FASTA',
                              'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                              'degen_translations': 'normal',
                              'exclude': 'YES',
                              'aminoacids': True,
                              'special': False,
                              'outgroup': '',
                          }
                          )
        html_page = res.content.decode('utf-8')
        expected = "Codon &#39;--C&#39; is invalid"
        self.assertFalse(expected in html_page)
