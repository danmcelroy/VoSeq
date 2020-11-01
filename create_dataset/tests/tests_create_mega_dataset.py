from django.test import TestCase
from django.conf import settings
from django.test.client import Client
from django.core.management import call_command

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreateMEGADatasetTest1(TestCase):
    def setUp(self):
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
            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
            'number_genes': None,
            'positions': ['ALL'],
            'degen_translations': None,
            'translations': False,
            'partition_by_positions': 'by gene',
            'file_format': 'MEGA',
            'aminoacids': False,
            'outgroup': '',
        }

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_headers(self):
        expected = '#MEGA\n!TITLE title;\n\n'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_sequence_line_breaks(self):
        expected = '#CP100_10_Aus_aus\nCGACGACGACGACGACGAC'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_sequence_concatenation(self):
        expected = 'CGACGACGACGACGACGACGAC'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_dataset_as_aminoacids(self):
        cleaned_data = self.cleaned_data
        cleaned_data['aminoacids'] = True
        dataset_creator = CreateDataset(cleaned_data)

        expected = 'DDDDDDDDDDDDDDDDDDDDDDDDDDD'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_dataset_with_partitions(self):
        cleaned_data = self.cleaned_data
        cleaned_data['partition_by_positions'] = '1st-2nd, 3rd'
        dataset_creator = CreateDataset(cleaned_data)

        expected = ''
        result = dataset_creator.dataset_str
        self.assertEqual(expected, result)

    def test_dataset_with_degen_translations(self):
        cleaned_data = self.cleaned_data
        cleaned_data['degen_translations'] = 'S'
        cleaned_data['translations'] = True
        dataset_creator = CreateDataset(cleaned_data)

        expected = 'GAYGAYGAYGAYGAYGAYGAYGAY'
        result = dataset_creator.dataset_str.strip()
        self.assertTrue(expected in result)

    def test_dataset_with_partitions_and_degen_translations(self):
        cleaned_data = self.cleaned_data
        cleaned_data['partition_by_positions'] = '1st-2nd, 3rd'
        cleaned_data['translations'] = True
        cleaned_data['degen_translations'] = 'normal'
        dataset_creator = CreateDataset(cleaned_data)

        expected = ''
        result = dataset_creator.dataset_str.strip()
        self.assertEqual(expected, result)
