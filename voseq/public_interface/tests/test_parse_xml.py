from django.test import TestCase

from public_interface.management.commands.migrate_db import ParseXML


class TestParseXML(TestCase):
    def setUp(self):
        with open('test_db_dump.xml') as handle:
            self.xml_string = handle.read()

        self.parse_xml_to_fail = ParseXML(xml_string=self.xml_string, tables_prefix='voseq_',
                                          verbosity=0)
        self.parse_xml = ParseXML(xml_string=self.xml_string, tables_prefix='',
                                  verbosity=0)

    def test_parse_table_genes(self):
        """We don't have tables with prefix 'voseq_' in our test_db_dump.xml file
        """
        self.assertRaises(ValueError, self.parse_xml_to_fail.parse_table_genes, self.xml_string)
