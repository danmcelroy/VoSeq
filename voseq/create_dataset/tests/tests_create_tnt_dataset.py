from django.test import TestCase
from django.test.client import Client
from django.core.management import call_command

from create_dataset.utils import CreateDataset
from public_interface.models import Genes


class CreateTNTDatasetTest(TestCase):
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
            'degen_translations': None,
            'translations': False,
            'positions': ['ALL'],
            'partition_by_positions': 'by gene',
            'file_format': 'TNT',
            'aminoacids': False,
            'outgroup': '',
        }

        self.c = Client()
        self.dataset_creator = CreateDataset(self.cleaned_data)
        self.maxDiff = None

    def test_create_dataset(self):
        expected = 'nstates dna;\nxread\n1909 2\n\n&[dna]\nCP100-10_Aus_aus'
        result = self.dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_create_dataset_outgroup(self):
        cleaned_data = self.cleaned_data
        cleaned_data['outgroup'] = 'CP100-11'
        dataset_creator = CreateDataset(cleaned_data)
        expected = 'nstates dna;\nxread\n1909 2\n\n&[dna]\nCP100-11_Aus_bus'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_create_dataset_all_codons(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['ALL']
        dataset_creator = CreateDataset(cleaned_data)
        expected = 'nstates dna;\nxread\n1909 2\n\n&[dna]\nCP100-10_Aus_aus'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

        expected = '??????????????????????????????????????????????????????????ATTATTCGAACAGAATTAAGTACCCCTGGATCATTAATCGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATGGTTATACCTATTATAATTGGAGGATTTGGTAATTGACTTATTCCCCTTATATTAGGAGCCCCTGATATAGCTTTTCCACGAATAAATAATATAAGATTTTGACTTCTCCCACCCTCTTTAATTTTATTAATTTCGAGTAGTATAGTAGAAAATGGTGCTGGCACAGGATGAACGGTCTATCCCCCCCTCTCATCTAATATTGCCCATAGAGGATCCTCAGTTGATTTAGCAATCTTTTCCTTACATTTAGCTGGAATCTCATCAATTCTTGGAGCAATTAATTTTATTACAACAATTATTAATATACGAATTAATAAAATATCTTATGATCAAATACCTTTATTTGTTTGAGCTGTAGGAATTACCGCATTATTATTATTACTTTCTTTACCTGTATTAGCTGGAGCTATCACAATACTACTCACAGATCGAAACTTAAATACATCTTTTTTTGACCCAGCAGGAGGTGGAGATCCTATTTTATATCAACATTTATTTTGATTTTTTGG'
        self.assertTrue(expected in result)

        expected = 'CP100-11_Aus_bus                                       ????????????????????????'
        self.assertTrue(expected in result)

    def test_create_dataset_1st_codon(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['1st']
        dataset_creator = CreateDataset(cleaned_data)
        expected = 'nstates dna;\nxread\n636 2\n\n&[dna]\nCP100-10_Aus_aus'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

        expected = '???????????????????AACAGTAACGTTAGGGCATAAAGAGCGTAAATTAGACAAAGGTGATCACCATG'
        self.assertTrue(expected in result)

    def test_fill_seqs_with_missing_chars(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['ALL']
        dataset_creator = CreateDataset(cleaned_data)
        result = dataset_creator.dataset_str
        expected = "CP100-10_Aus_aus                                       ??????????????????????????????????????????????????????????ATTATTCGAACAGAATTAAGTAC"
        self.assertTemplateNotUsed(expected in result)

    def test_create_dataset_aa(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['ALL']
        cleaned_data['aminoacids'] = True
        dataset_creator = CreateDataset(cleaned_data)
        expected = 'XXXXXXXXXXXXXXXXXXXIIRTELSTPGSLI'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

        expected = 'xread\n635 2\n'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_create_dataset_aa_with_outgroup(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['ALL']
        cleaned_data['outgroup'] = 'CP100-11'
        cleaned_data['aminoacids'] = True
        dataset_creator = CreateDataset(cleaned_data)
        expected = '&[protein]\nCP100-11_Aus_bus                                       TLYFIFGIWAGM'
        result = dataset_creator.dataset_str
        self.assertTrue(expected in result)

    def test_create_dataset_try_with_degen_translation_and_paritioned_by_codon(self):
        cleaned_data = self.cleaned_data
        cleaned_data['positions'] = ['ALL']
        cleaned_data['degen_translations'] = 'normal'
        cleaned_data['translations'] = True
        cleaned_data['partition_by_positions'] = '1st-2nd, 3rd'
        dataset_creator = CreateDataset(cleaned_data)
        expected = ''

        result = dataset_creator.dataset_str
        self.assertEqual(expected, result)
