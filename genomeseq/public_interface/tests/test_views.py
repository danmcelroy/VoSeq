import json
import os

from django.test import Client
from django.test import TestCase

from public_interface.models import Vouchers


class TestViews(TestCase):
    def setUp(self):
        json_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "NW1-1.json")
        with open(json_file, "r") as handle:
            item = json.loads(handle.read())
            item['max_altitude'] = None
            item['min_altitude'] = None
        b = Vouchers.objects.create(**item)
        b.save()

    def test_index(self):
        c = Client()
        response = c.get('/')
        self.assertEqual(200, response.status_code)

    def test_browse(self):
        c = Client()
        response = c.get('/browse')
        self.assertEqual(200, response.status_code)

    def test_show_voucher(self):
        c = Client()
        response = c.get('/p/NW1-1', follow=True)
        self.assertEqual(200, response.status_code)
