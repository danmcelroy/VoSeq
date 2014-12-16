import json
import os

from django.test import Client
from django.test import TestCase

from public_interface.models import Vouchers
from public_interface.models import FlickrImages
from public_interface.models import Sequences


class TestViews(TestCase):
    def setUp(self):
        json_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "NW1-1.json")
        json_file_seqs = os.path.join(os.path.dirname(os.path.abspath(__file__)), "NW1-1_seqs.json")
        with open(json_file_seqs, "r") as handle:
            seqs = json.loads(handle.read())

        with open(json_file, "r") as handle:
            item = json.loads(handle.read())
            item['max_altitude'] = None
            item['min_altitude'] = None
        b = Vouchers.objects.create(**item)
        b.save()

        f = FlickrImages.objects.create(voucher=b)
        f.save()

        s = Sequences.objects.create(code=b, sequences=seqs['sequences'])
        s.save()

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

    def test_show_voucher_doesnt_exist(self):
        c = Client()
        response = c.get('/p/NW1-1aaaaaaa', follow=True)
        self.assertEqual(404, response.status_code)
