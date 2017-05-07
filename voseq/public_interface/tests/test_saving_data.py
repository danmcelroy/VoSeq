from django.core.management import call_command
from django.db import connection
from django.test import TestCase

from public_interface.models import Sequences
from public_interface.models import Vouchers


class TestViews(TestCase):
    def setUp(self):
        with connection.cursor() as cursor:
            cursor.execute("alter sequence public_interface_vouchers_id_seq restart with 1")
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        self.voucher_model = Vouchers.objects.filter(code='CP100-18')[0]

    def test_save_sequences_ambiguous_characters(self):
        sequence_model = Sequences(
            code=self.voucher_model,
            sequences='???---NNNATCTACTA',
            gene_code='COI',
            genbank=False,
        )
        sequence_model.save()

        sequence_model = Sequences.objects.filter(
            code=self.voucher_model,
            gene_code='COI',
        )[0]
        self.assertEqual(sequence_model.number_ambiguous_bp, 25)
