import os

from django.conf import settings
from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command

from create_dataset.utils import CreateDataset
from public_interface.models import Genes
from public_interface.models import GeneSets
from public_interface.models import TaxonSets


class CreateMegaDatasetTest(TestCase):
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
            'degen_translations': None,
            'number_genes': None,
            'positions': ['ALL'],
            'translations': False,
            'partition_by_positions': 'by gene',
            'file_format': 'MEGA',
            'aminoacids': False,
            'outgroup': '',
        }

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_all_codons_as_one(self):
        dataset_file = os.path.join(settings.BASE_DIR, '..', 'create_dataset',
                                    'tests', 'create_mega_dataset', 'dataset.meg')
        with open(dataset_file, 'r') as handle:
            expected = handle.read()

        cleaned_data = self.cleaned_data.copy()
        cleaned_data['positions'] = ['ALL']
        cleaned_data['partition_by_positions'] = 'by gene'
        dataset_creator = CreateDataset(cleaned_data)

        result = dataset_creator.dataset_str
        self.assertEqual(expected.strip(), result)

    def test_all_codons_partitioned_as_each(self):
        cleaned_data = self.cleaned_data.copy()
        cleaned_data['positions'] = ['ALL']
        cleaned_data['partition_by_positions'] = 'by codon position'
        dataset_creator = CreateDataset(cleaned_data)
        expected = 'Cannot produce MEGA dataset with codon positions in different partitions.'
        result = dataset_creator.errors
        self.assertTrue(expected in result)
