from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command
from django.contrib.auth.models import User

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreateFASTADatasetTest(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump2.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        g1 = Genes.objects.get(gene_code='COI-begin')
        g2 = Genes.objects.get(gene_code='ef1a')
        self.cleaned_data = {
            'gene_codes': [g1, g2],
            'taxonset': None,
            'voucher_codes': 'CP100-10\r\nCP100-11',
            'geneset': None,
            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
            'number_genes': None,
            'translations': False,
            'degen_translations': 'normal',
            'positions': ['ALL'],
            'partition_by_positions': 'by gene',
            'file_format': 'FASTA',
            'aminoacids': True,
            'outgroup': '',
        }

        self.user = User.objects.get(username='admin')
        self.user.set_password('pass')
        self.user.save()

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_create_dataset_degenerated(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 3,  # wingless
                            'geneset': '',
                            'taxonset': '',
                            'translations': True,
                            'introns': 'YES',
                            'file_format': 'FASTA',
                            'degen_translations': 'normal',
                            'exclude': 'YES',
                            'aminoacids': False,
                            'special': False,
                            'outgroup': '',
                            'positions': 'ALL',
                            'partition_by_positions': 'by gene',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = 'TNGGNTTYATHGTNTGAGCNCAYCAYATHTTYACN'
        self.assertTrue(expected in str(c.content))

    def test_create_dataset_degenerated_warning_data_cannot_be_partitioned(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 4,
                            'geneset': '',
                            'taxonset': '',
                            'introns': 'YES',
                            'file_format': 'FASTA',
                            'translations': True,
                            'degen_translations': 'normal',
                            'exclude': 'YES',
                            'aminoacids': False,
                            'special': False,
                            'outgroup': '',
                            'positions': 'ALL',
                            'partition_by_positions': 'by codon position',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = 'Cannot degenerate codons if they go to different partitions'
        self.assertTrue(expected in str(c.content))

    def test_create_dataset_degenerated_warning_data_cannot_be_of_partial_codons(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 4,
                            'geneset': '',
                            'taxonset': '',
                            'introns': 'YES',
                            'file_format': 'FASTA',
                            'translations': True,
                            'degen_translations': 'normal',
                            'exclude': 'YES',
                            'aminoacids': False,
                            'special': False,
                            'outgroup': '',
                            'positions': '1st',
                            'partition_by_positions': 'by gene',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = 'Cannot degenerate codons if you have not selected all codon positions'
        self.assertTrue(expected in str(c.content))

    def test_fasta_as_aminoacids(self):
        self.c.post('/accounts/login/', {'username': 'admin', 'password': 'pass'})
        c = self.c.post('/create_dataset/results/',
                        {
                            'voucher_codes': 'CP100-10',
                            'gene_codes': 3,  # wingless
                            'geneset': '',
                            'taxonset': '',
                            'translations': True,
                            'introns': 'YES',
                            'file_format': 'FASTA',
                            'degen_translations': 'normal',
                            'exclude': 'YES',
                            'aminoacids': True,
                            'special': False,
                            'outgroup': '',
                            'positions': 'ALL',
                            'partition_by_positions': 'by gene',
                            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
                        }
                        )
        expected = 'IYAMLAIGLLGFIVWAHHM'
        self.assertTrue(expected in str(c.content))
