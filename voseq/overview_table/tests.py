from django.test import TestCase
from django.core.management import call_command

from overview_table.utils import OverviewTableMaker


class TestOverviewTable(TestCase):

    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_utils(self):
        o = OverviewTableMaker()
        expected = 3
        result = len(o.items)
        self.assertEqual(expected, result)
