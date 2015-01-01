from django.core.management import call_command
from django.test import TestCase

from public_interface.models import Vouchers


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
