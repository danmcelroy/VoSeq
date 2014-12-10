#!-*- encoding: utf-8 -*-
"""
Needs an XML file with the database dump from MySQL:

> mysqldump --xml database > dump.xml
"""
import codecs
import datetime
import dataset
import json
import re
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
        self.table_genes_items = None
        self.table_genesets_items = None
        self.table_members_items = None
        self.table_primers_items = None
        self.table_sequences_items = None
        self.table_taxonsets_items = None
        self.table_vouchers_items = None

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

    def parse_table_sequences(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "sequences"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_sequences_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['code'] = row.find("./field/[@name='code']").text
            item['geneCode'] = row.find("./field/[@name='geneCode']").text
            item['sequences'] = row.find("./field/[@name='sequences']").text
            item['accession'] = row.find("./field/[@name='accession']").text
            item['labPerson'] = row.find("./field/[@name='labPerson']").text
            item['dateCreation'] = row.find("./field/[@name='dateCreation']").text
            item['dateModification'] = row.find("./field/[@name='dateModification']").text
            item['notes'] = row.find("./field/[@name='notes']").text
            item['genbank'] = row.find("./field/[@name='genbank']").text
            item['timestamp'] = row.find("./field/[@name='timestamp']").text
            self.table_sequences_items.append(item)

    def parse_table_taxonsets(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "taxonsets"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_taxonsets_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['taxonset_name'] = row.find("./field/[@name='taxonset_name']").text
            item['taxonset_creator'] = row.find("./field/[@name='taxonset_creator']").text
            item['taxonset_description'] = row.find("./field/[@name='taxonset_description']").text
            item['taxonset_list'] = row.find("./field/[@name='taxonset_list']").text
            item['taxonset_id'] = row.find("./field/[@name='taxonset_id']").text
            self.table_taxonsets_items.append(item)

    def parse_table_vouchers(self, xml_string):
        our_data = False
        this_table = self.tables_prefix + "vouchers"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        if our_data is False:
            raise ValueError("Could not find table %s in database dump file." % this_table)

        self.table_vouchers_items = []
        for row in our_data.findall('row'):
            item = dict()
            item['code'] = row.find("./field/[@name='code']").text
            item['orden'] = row.find("./field/[@name='orden']").text
            item['family'] = row.find("./field/[@name='family']").text
            item['subfamily'] = row.find("./field/[@name='subfamily']").text
            item['tribe'] = row.find("./field/[@name='tribe']").text
            item['subtribe'] = row.find("./field/[@name='subtribe']").text
            item['genus'] = row.find("./field/[@name='genus']").text
            item['species'] = row.find("./field/[@name='species']").text
            item['subspecies'] = row.find("./field/[@name='subspecies']").text
            item['country'] = row.find("./field/[@name='country']").text
            item['specificLocality'] = row.find("./field/[@name='specificLocality']").text
            item['typeSpecies'] = row.find("./field/[@name='typeSpecies']").text
            item['latitude'] = row.find("./field/[@name='latitude']").text
            item['longitude'] = row.find("./field/[@name='longitude']").text
            item['altitude'] = row.find("./field/[@name='altitude']").text
            item['collector'] = row.find("./field/[@name='collector']").text
            item['dateCollection'] = row.find("./field/[@name='dateCollection']").text
            item['voucherImage'] = row.find("./field/[@name='voucherImage']").text
            item['thumbnail'] = row.find("./field/[@name='thumbnail']").text
            item['extraction'] = row.find("./field/[@name='extraction']").text
            item['dateExtraction'] = row.find("./field/[@name='dateExtraction']").text
            item['extractor'] = row.find("./field/[@name='extractor']").text
            item['voucherLocality'] = row.find("./field/[@name='voucherLocality']").text
            item['publishedIn'] = row.find("./field/[@name='publishedIn']").text
            item['notes'] = row.find("./field/[@name='notes']").text
            item['edits'] = row.find("./field/[@name='edits']").text
            item['latesteditor'] = row.find("./field/[@name='latesteditor']").text
            item['hostorg'] = row.find("./field/[@name='hostorg']").text
            item['sex'] = row.find("./field/[@name='sex']").text
            item['extractionTube'] = row.find("./field/[@name='extractionTube']").text
            item['voucher'] = row.find("./field/[@name='voucher']").text
            item['voucherCode'] = row.find("./field/[@name='voucherCode']").text
            item['flickr_id'] = row.find("./field/[@name='flickr_id']").text
            item['determinedBy'] = row.find("./field/[@name='determinedBy']").text
            item['auctor'] = row.find("./field/[@name='auctor']").text
            item['timestamp'] = row.find("./field/[@name='timestamp']").text
            self.table_vouchers_items.append(item)

    def convert_to_int(self, string):
        try:
            string = int(string)
        except TypeError:
            string = None
        except ValueError:
            string = None
        return string

    def import_table_vouchers(self):
        if self.table_vouchers_items is None:
            self.parse_table_vouchers(self.dump_string)

        for item in self.table_vouchers_items:
            if item['altitude'] is not None:
                altitude = re.sub("\s+", "", item['altitude'])
                altitude = altitude.split("-")

                if len(altitude) > 1:
                    max_altitude = altitude[0].strip()
                    max_altitude = re.sub("[a-zA-Z]", "", max_altitude)
                    max_altitude = self.convert_to_int(max_altitude)

                    min_altitude = altitude[1].strip()
                    min_altitude = re.sub("[a-zA-Z]", "", min_altitude)
                    min_altitude = self.convert_to_int(min_altitude)

                    item['max_altitude'] = max_altitude
                    item['min_altitude'] = min_altitude
                else:
                    max_altitude = re.sub("[a-zA-Z]", "", altitude[0])
                    max_altitude = re.sub("[a-zA-Z]", "", max_altitude)
                    max_altitude = self.convert_to_int(max_altitude)
                    item['max_altitude'] = max_altitude
                    item['min_altitude'] = None
            else:
                item['max_altitude'] = None
                item['min_altitude'] = None

            print(item['max_altitude'], item['min_altitude'])

dump_file = sys.argv[1].strip()
with codecs.open(dump_file, "r") as handle:
    dump = handle.read()

# tables_prefix = 'voseq_'
tables_prefix = ''
parser = ParseXML(dump, tables_prefix)
#print(parser.table_genes_items)
#print(parser.table_genesets_items)
#print(parser.table_members_items)
#print(parser.table_primers_items)
#print(parser.table_sequences_items)
# print(parser.table_taxonsets_items)
parser.import_table_vouchers()
