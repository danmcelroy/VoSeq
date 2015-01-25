from django.test import TestCase
from django.core.management import call_command

from stats.models import Stats


class TestCustomCommand(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_create_stats(self):
        call_command('create_stats')
        res = Stats.objects.get(id=1)
        expected = 10
        self.assertEqual(expected, res.vouchers)
