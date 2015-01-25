from django.core.management import call_command
from django.test import Client
from django.test import TestCase


class TestViews(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        self.client = Client()

    def test_genes(self):
        response = self.client.get('/genes/')
        self.assertEqual(200, response.status_code)

    def test_view_gene(self):
        response = self.client.get('/genes/wingless/')
        self.assertEqual(200, response.status_code)

    def test_view_gene_not_found(self):
        response = self.client.get('/genes/winglessaaaaaa/')
        self.assertEqual(200, response.status_code)
