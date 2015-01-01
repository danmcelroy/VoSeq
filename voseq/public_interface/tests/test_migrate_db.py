from django.core.management import call_command
from django.test import TestCase

from public_interface.models import Vouchers


class TestCustomCommand(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml'}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_superfamily(self):
        b = Vouchers.objects.get(code='CP100-10')
        self.assertEqual('', b.superfamily)
