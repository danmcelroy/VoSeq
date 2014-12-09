#!-*- encoding: utf-8 -*-
"""
Needs an XML file with the database dump from MySQL:

> mysqldump --xml database > dump.xml
"""
import codecs
import datetime
import dataset
import json
import sys
import xml.etree.ElementTree as ET


if len(sys.argv) < 2:
    print("Enter name of database dump file as argument.")
    print("This file can be obtained from your MySQL database using this command")
    print("\t> mysqdump --xml database > dump.xml")
    sys.exit(1)

with open("config.json", "r") as f:
    settings = json.loads(f.read())


db_url = 'postgresql://' + settings['DB_USER'] + ':' + settings['DB_PASS'] + '@' \
             + settings['DB_HOST'] + ":" + settings['DB_PORT'] + "/" + settings['DB_NAME']
db = dataset.connect(db_url)


class ParseXML(object):
    """
    Parses MySQL dump as XML file.
    """
    def __init__(self, xml_string, tables_prefix=None):
        if tables_prefix is None:
            self.tables_prefix = ''
        else:
            self.tables_prefix = tables_prefix

        self.dump_string = xml_string
        self.table_genes = self.parse_table_genes(xml_string)
        self.table_genesets = self.parse_table_genesets(xml_string)
        self.table_members = self.parse_table_members(xml_string)
        self.table_primers = self.parse_table_primers(xml_string)

    def parse_table_genes(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "genes"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_genes_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['geneCode'] = row.find("./field/[@name='geneCode']").text
            item['length'] = row.find("./field/[@name='length']").text
            item['description'] = row.find("./field/[@name='description']").text
            item['readingframe'] = row.find("./field/[@name='readingframe']").text
            item['notes'] = row.find("./field/[@name='notes']").text
            item['timestamp'] = row.find("./field/[@name='timestamp']").text
            item['genetic_code'] = row.find("./field/[@name='genetic_code']").text
            item['aligned'] = row.find("./field/[@name='aligned']").text
            item['intron'] = row.find("./field/[@name='intron']").text
            item['prot_code'] = row.find("./field/[@name='prot_code']").text
            item['genetype'] = row.find("./field/[@name='genetype']").text
            self.table_genes_items.append(item)

    def parse_table_genesets(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "genesets"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_genesets_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['geneset_name'] = row.find("./field/[@name='geneset_name']").text
            item['geneset_creator'] = row.find("./field/[@name='geneset_creator']").text
            item['geneset_description'] = row.find("./field/[@name='geneset_description']").text
            item['geneset_list'] = row.find("./field/[@name='geneset_list']").text
            item['geneset_id'] = row.find("./field/[@name='geneset_id']").text
            self.table_genesets_items.append(item)

    def parse_table_members(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "members"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_members_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['member_id'] = row.find("./field/[@name='member_id']").text
            item['firstname'] = row.find("./field/[@name='firstname']").text
            item['lastname'] = row.find("./field/[@name='lastname']").text
            item['login'] = row.find("./field/[@name='login']").text
            item['passwd'] = row.find("./field/[@name='passwd']").text
            item['admin'] = row.find("./field/[@name='admin']").text
            self.table_members_items.append(item)

    def parse_table_primers(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "primers"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_primers_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['code'] = row.find("./field/[@name='code']").text
            item['geneCode'] = row.find("./field/[@name='geneCode']").text
            item['primer1'] = row.find("./field/[@name='primer1']").text
            item['primer2'] = row.find("./field/[@name='primer2']").text
            item['primer3'] = row.find("./field/[@name='primer3']").text
            item['primer4'] = row.find("./field/[@name='primer4']").text
            item['primer5'] = row.find("./field/[@name='primer5']").text
            item['primer6'] = row.find("./field/[@name='primer6']").text
            item['timestamp'] = row.find("./field/[@name='timestamp']").text
            self.table_primers_items.append(item)


dump_file = sys.argv[1].strip()
with codecs.open(dump_file, "r") as handle:
    dump = handle.read()

# tables_prefix = 'voseq_'
tables_prefix = ''
parser = ParseXML(dump, tables_prefix)
#print(parser.table_genes_items)
#print(parser.table_genesets_items)
#print(parser.table_members_items)
print(parser.table_primers_items)
"""
print(parser.table_sequences_items)
print(parser.table_taxonsets_items)
print(parser.table_vouchers_items)
"""
