from django.test import TestCase
from django.core.management import call_command

from public_interface.models import GeneSets
from public_interface.models import TaxonSets
from voucher_table.views import VoucherTable


class TestVoucherTable(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        taxonset = TaxonSets.objects.all()[0]
        geneset = GeneSets.objects.all()[0]

        self.cleaned_data = {
            'collector_info': ['country', 'specificLocality', 'collector'],
            'gene_codes': [],
            'gene_info': 'NUMBER OF BASES',
            'voucher_info': ['code', 'genus', 'species'],
            'field_delimitor': 'COMMA',
            'taxonset': taxonset,
            'voucher_codes': '',
            'geneset': geneset,
        }
        self.table = VoucherTable(self.cleaned_data)

    def test_header_csv_file(self):
        expected = ['Code', 'Genus', 'Species', 'Country', 'Specific Locality', 'Collector', 'wingless',
                    '16S', 'EF1a', 'COI']
        result = self.table.get_headers()
        self.assertEqual(expected, result)
