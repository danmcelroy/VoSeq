import json
import os

from django.core.management import call_command
from django.test import Client
from django.test import TestCase
from django.test.utils import override_settings
import haystack

from public_interface.models import Vouchers
from public_interface.models import FlickrImages
from public_interface.models import Sequences
from public_interface import views


# Need to use a clean index for our tests
TEST_INDEX = {
    'default': {
        'ENGINE': 'haystack.backends.elasticsearch_backend.ElasticsearchSearchEngine',
        'URL': 'http://127.0.0.1:9200/',
        'TIMEOUT': 60 * 10,
        'INDEX_NAME': 'test_index',
    },
}


@override_settings(HAYSTACK_CONNECTIONS=TEST_INDEX)
class TestViews(TestCase):
    def setUp(self):
        json_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "public_interface", "vouchers.json")
        json_file_seqs = os.path.join(os.path.dirname(os.path.abspath(__file__)), "public_interface", "NN1-1_seqs.json")
        with open(json_file_seqs, "r") as handle:
            seqs = json.loads(handle.read())

        with open(json_file, "r") as handle:
            items = json.loads(handle.read())

            for item in items:
                item['max_altitude'] = None
                item['min_altitude'] = None
                b = Vouchers.objects.create(**item)
                b.save()

                f = FlickrImages.objects.create(voucher=b)
                f.save()

                s = Sequences.objects.create(code=b, sequences=seqs['sequences'])
                s.save()

        # build index with our test data
        haystack.connections.reload('default')
        call_command('rebuild_index', interactive=False, verbosity=0)
        super(TestViews, self).setUp()

        self.client = Client()

    def test_index(self):
        response = self.client.get('/')
        self.assertEqual(200, response.status_code)

    def test_browse(self):
        response = self.client.get('/browse/')
        self.assertEqual(200, response.status_code)

    def test_show_voucher(self):
        response = self.client.get('/p/NN1-1', follow=True)
        self.assertEqual(200, response.status_code)

    def test_show_voucher_doesnt_exist(self):
        response = self.client.get('/p/NN1-1aaaaaaa', follow=True)
        self.assertEqual(404, response.status_code)

    def test_search_hymenoptera(self):
        response = self.client.get('/search/?orden=Hymenoptera')
        content = str(response.content)
        if 'NN1-2' in content and 'NN1-1' not in content:
            found_item = True
        else:
            found_item = False
        self.assertTrue(found_item)

    def tearDown(self):
        call_command('clear_index', interactive=False, verbosity=0)
