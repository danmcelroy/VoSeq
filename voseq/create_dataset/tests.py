from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreateDatasetTest(TestCase):
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
        }

    def test_create_dataset(self):
        dataset_creator = CreateDataset(self.cleaned_data)
        expected = '>CP100-11\n??TGAGCCGGTATAATTGGTACATCCCTAAGTCTTATTATTC'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_create_dataset_invalid_voucher(self):
        self.cleaned_data['voucher_codes'] = 'CP100-10\r\nCP100-11\r\nCP555-55'

        expected = 'Could not find voucher CP555-55'
        dataset_creator = CreateDataset(self.cleaned_data)
        result = dataset_creator.errors
        self.assertTrue(expected in result)
