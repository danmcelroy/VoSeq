from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreateMEGADatasetTest(TestCase):
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
            'positions': ['ALL'],
            'degen_translations': None,
            'partition_by_positions': 'ONE',
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
        expected = '#CP100-10_Melitaea_diamina\n?????????????????????????TGAGCCGGTATAATTGGTACAT'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_sequence_concatenation(self):
        expected = '??????????????????????????CAAGT'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_dataset_as_aminoacids(self):
        cleaned_data = self.cleaned_data
        cleaned_data['aminoacids'] = True
        dataset_creator = CreateDataset(cleaned_data)

        expected = 'PSFLIGDDQIYNTIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLG'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_dataset_with_partitions(self):
        cleaned_data = self.cleaned_data
        cleaned_data['partition_by_positions'] = '1st2nd_3rd'
        dataset_creator = CreateDataset(cleaned_data)

        expected = '?TGGCGGATATGGACTCCTAGCTATATCGACGATTGGAACCAGTTTTATGGGAGAC'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_dataset_with_degen_tranlations(self):
        cleaned_data = self.cleaned_data
        cleaned_data['degen_translations'] = 'NORMAL'
        cleaned_data['translations'] = True
        dataset_creator = CreateDataset(cleaned_data)

        expected = 'NTGRGCNGGNATRATYGGNA'
        result = dataset_creator.dataset_str.strip()
        self.assertTrue(expected in result)

    def test_dataset_with_partitions_and_degen_tranlations(self):
        cleaned_data = self.cleaned_data
        cleaned_data['partition_by_positions'] = '1st2nd_3rd'
        cleaned_data['translations'] = True
        dataset_creator = CreateDataset(cleaned_data)

        expected = '#MEGA\n!TITLE title;'
        result = dataset_creator.dataset_str.strip()
        self.assertEqual(expected, result)
