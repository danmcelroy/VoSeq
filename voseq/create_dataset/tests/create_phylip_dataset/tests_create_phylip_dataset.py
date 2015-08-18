from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command
from django.contrib.auth.models import User

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreatePhylipDatasetTest(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        g1 = Genes.objects.get(gene_code='COI')
        g2 = Genes.objects.get(gene_code='EF1a')
        self.cleaned_data = {
            'gene_codes': [g1, g2],
            'taxonset': None,
            'voucher_codes': 'CP100-10\r\nCP100-11',
            'geneset': None,
            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
            'number_genes': None,
            'degen_translations': None,
            'positions': ['ALL'],
            'partition_by_positions': 'ONE',
            'file_format': 'PHY',
            'aminoacids': True,
            'outgroup': '',
        }

        self.user = User.objects.get(username='admin')
        self.user.set_password('pass')
        self.user.save()

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_numer_of_chars_for_aa_dataset(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 3,  # wingless
                            'geneset': '',
                            'taxonset': '',
                            'translations': False,
                            'introns': 'YES',
                            'file_format': 'PHY',
                            'degen_translations': 'NORMAL',
                            'exclude': 'YES',
                            'aminoacids': True,
                            'special': False,
                            'outgroup': '',
                            'positions': 'ALL',
                            'partition_by_positions': 'ONE',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = '1 137'
        self.assertTrue(expected in str(c.content))

    def test_partitioned_1st2nd_3rd(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 3,  # wingless
                            'geneset': '',
                            'taxonset': '',
                            'translations': False,
                            'introns': 'YES',
                            'file_format': 'PHY',
                            'degen_translations': 'NORMAL',
                            'exclude': 'YES',
                            'aminoacids': False,
                            'special': False,
                            'outgroup': '',
                            'positions': 'ALL',
                            'partition_by_positions': '1st2nd_3rd',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = 'CAGTGATCGGAATCACACACGGCATTATT'
        self.assertTrue(expected in str(c.content))
