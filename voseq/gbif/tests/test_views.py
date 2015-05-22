from django.test import TestCase
from django.test import Client


class TestViews(TestCase):
    def setUp(self):
        self.client = Client()

    def test_index(self):
        response = self.client.get('/share_data_gbif/')
        self.assertEqual(200, response.status_code)

    def test_dump_data_error(self):
        response = self.client.get('/share_data_gbif/dump_data/')
        result = response.content.decode('utf-8')
        expected = '{"result": "error"}'
        self.assertEqual(expected, result)
