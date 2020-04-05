from django.test import TestCase
from django.core.management import call_command

from core import exceptions
from core.utils import clean_positions
from core.utils import get_gene_codes
from core.utils import get_voucher_codes
from public_interface.models import TaxonSets
from public_interface.models import GeneSets
from public_interface.models import Genes


class TestCoreUtils(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        gs = GeneSets.objects.get(geneset_name='2genes')
        g = Genes.objects.get(gene_code='COI')
        g2 = Genes.objects.get(gene_code='16S')
        ts = TaxonSets.objects.get(taxonset_name='Erebia')
        self.cleaned_data = {
            'gene_codes': [g, g2],
            'taxonset': ts,
            'voucher_codes': 'CP200-10\r\n \r\nCP100-11\r\n \r\n',
            'geneset': gs,
        }

    def test_get_gene_codes(self):
        expected = 3
        result = get_gene_codes(self.cleaned_data)
        self.assertEqual(expected, len(result))

    def test_get_voucher_codes(self):
        expected = 3
        result = get_voucher_codes(self.cleaned_data)
        self.assertEqual(expected, len(result))

    def test_get_voucher_codes_order(self):
        """Voucher codes should not be sorted. Same order as input from author
        should be kept."""
        self.cleaned_data['voucher_codes'] = 'CP100-10\r\nCP100-11\r\nCP100-12\r\nCP100-13\r\nCP100-11'
        expected = ('CP100-10', 'CP100-11', 'CP100-12', 'CP100-13')
        result = get_voucher_codes(self.cleaned_data)
        self.assertEqual(expected, result)

    def test_get_voucher_codes_dropped(self):
        self.cleaned_data['voucher_codes'] = 'CP100-10\r\n--CP100-11\r\nCP100-12'
        expected = 2
        result = get_voucher_codes(self.cleaned_data)
        self.assertEqual(expected, len(result))

    def test_clean_positions(self):
        self.assertEqual(clean_positions(['ALL', '1st']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '2nd']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '3rd']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '1st', '2nd']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '1st', '3rd']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '2nd', '3rd']), ['ALL'], 'Has "ALL" in list.')
        self.assertEqual(clean_positions(['ALL', '1st', '2nd', '3rd']), ['ALL'], 'Has "ALL" in list.')

        self.assertEqual(clean_positions(['1st', '2nd', '3rd']), ['ALL'], 'All codon positions were required.')
        self.assertEqual(clean_positions(['1st', '2nd']), ['1st-2nd'], 'Only some codon positions were required.')
        self.assertRaises(exceptions.InadequateCodonPositions, clean_positions, ['1st', '3rd'])
        self.assertRaises(exceptions.InadequateCodonPositions, clean_positions, ['2nd', '3rd'])

        self.assertEqual(clean_positions(['1st']), ['1st'], 'All codon positions were required.')
        self.assertEqual(clean_positions(['2nd']), ['2nd'], 'All codon positions were required.')
        self.assertEqual(clean_positions(['3rd']), ['3rd'], 'All codon positions were required.')
