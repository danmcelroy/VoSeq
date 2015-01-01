from django.core.management import call_command
from django.test import TestCase

from public_interface.models import Vouchers
from public_interface.models import Sequences


class TestCustomCommand(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml'}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_orden_none(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual('', b.orden)

    def test_orden_null(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual('', b.orden)

    def test_orden_space(self):
        b = Vouchers.objects.get(code='CP100-12')
        self.assertEqual('', b.orden)

    def test_superfamily(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual('', b.superfamily)

    def test_family(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual('Nymphalidae', b.family)

    def test_family_empty(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual('', b.family)

    def test_family_space(self):
        b = Vouchers.objects.get(code='CP100-11')
        self.assertEqual('', b.family)

    def test_subfamily(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual('Nymphalinae', b.subfamily)

    def test_subfamily_empty(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual('', b.subfamily)

    def test_subfamily_space(self):
        b = Vouchers.objects.get(code='CP100-11')
        self.assertEqual('', b.subfamily)

    def test_subfamily_null(self):
        b = Vouchers.objects.get(code='CP100-12')
        self.assertEqual('', b.subfamily)

    def test_tribe_null_with_space(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual('', b.tribe)

    def test_genus(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual('Melitaea?', b.genus)

    def test_accession(self):
        b = Sequences.objects.get(code='CP100-10', gene_code='EF1a')
        self.assertEqual('AY218269', b.accession)

    def test_max_altitude_null_as_str(self):
        b = Vouchers.objects.get(code='CP100-09')
        self.assertEqual(None, b.max_altitude)

    def test_max_altitude_null(self):
        b = Vouchers.objects.get(code='CP100-11')
        self.assertEqual(None, b.max_altitude)

    def test_max_altitude_empty_with_space(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual(None, b.max_altitude)

    def test_max_altitude1(self):
        b = Vouchers.objects.get(code='CP100-12')
        self.assertEqual(600, b.max_altitude)

    def test_max_altitude2(self):
        b = Vouchers.objects.get(code='CP100-13')
        self.assertEqual(2500, b.max_altitude)

    def test_max_altitude3(self):
        b = Vouchers.objects.get(code='CP100-14')
        self.assertEqual(2000, b.max_altitude)
