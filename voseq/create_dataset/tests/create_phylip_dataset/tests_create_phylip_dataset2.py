import os

from django.test import TestCase
from django.test.client import Client
from django.conf import settings
from django.core.management import call_command
from django.contrib.auth.models import User

from create_dataset.utils import CreateDataset
from public_interface.models import GeneSets
from public_interface.models import TaxonSets


class CreatePhylipDatasetTest(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump2.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        gene_set = GeneSets.objects.get(geneset_name='all_genes')
        taxon_set = TaxonSets.objects.get(taxonset_name='all_taxa')
        self.cleaned_data = {
            'gene_codes': '',
            'taxonset': taxon_set,
            'voucher_codes': '',
            'geneset': gene_set,
            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
            'number_genes': None,
            'degen_translations': None,
            'positions': ['ALL'],
            'partition_by_positions': 'ONE',
            'file_format': 'PHY',
            'aminoacids': False,
            'outgroup': '',
        }

        self.user = User.objects.get(username='admin')
        self.user.set_password('pass')
        self.user.save()

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_create_simple_dataset(self):
        dataset_file = os.path.join(settings.BASE_DIR, '..', 'create_dataset',
                                    'tests', 'create_phylip_dataset', 'dataset.phy')
        with open(dataset_file, "r") as handle:
            expected = handle.read()
        result = self.dataset_creator.dataset_str
        self.assertEqual(expected, result)

    def test_charset_block_file_of_simple_dataset(self):
        dataset_file = os.path.join(settings.BASE_DIR, '..', 'create_dataset',
                                    'tests', 'create_phylip_dataset', 'phylip_file.txt')
        with open(dataset_file, "r") as handle:
            expected = handle.read()
        result = self.dataset_creator.charset_block
        self.assertEqual(expected, result)
