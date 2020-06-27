import unittest

from public_interface.management.commands._migrate_db import get_as_tuple


class Test_Migrate(unittest.TestCase):
    def test_get_as_tuple(self):
        my_string = "|14799231204"
        result = get_as_tuple(my_string)
        expected = ("14799231204",)
        self.assertEqual(expected, result)
