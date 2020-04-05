from django.test import TestCase
from django.test import Client
from django.core.management import call_command


class TestViews(TestCase):
    def setUp(self):
        self.client = Client()

        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_index(self):
        response = self.client.get('/blast_local/CP100-10/COI/')
        self.assertEqual(200, response.status_code)
