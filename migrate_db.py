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


def migrate_genes_table(old_db, new_db):
    fixed_items = []
    append = fixed_items.append

    res = old_db.query("select * from genes")
    for i in res:
        del i['id']
        del i['timestamp']
        i['gene_code'] = i['geneCode']
        del i['geneCode']
        if i['genetic_code'] is None:
            i['genetic_code'] = 999
        i['reading_frame'] = i['readingframe']
        del i['readingframe']
        if i['reading_frame'] is None:
            i['reading_frame'] = 999

        if i['notes'] is None:
            i['notes'] = ''
        if i['intron'] is None:
            i['intron'] = ''

        i['gene_type'] = i['genetype']
        del i['genetype']
        if i['gene_type'] is None:
            i['gene_type'] = ''

        i['time_created'] = datetime.date.today()

        append(i)

    table = new_db['public_interface_genes']
    table.insert_many(fixed_items)


def migrate_genesets_table(old_db, new_db):
    print(old_db['genesets'].columns)
    fixed_items = []
    append = fixed_items.append
    res = old_db.query("select * from genesets")
    for i in res:
        del i['id']
        append(i)

    table = new_db['public_interface_genesets']
    table.insert_many(fixed_items)


def migrate_members_table(old_db, new_db):
    fixed_items = []
    append = fixed_items.append
    res = old_db.query("select * from members")
    for i in res:
        del i['id']
        i['firstname'] = i['firstname'].decode('utf-8')
        try:
            i['lastname'] = i['lastname'].decode('utf-8')
        except UnicodeDecodeError:
            i['lastname'] = i['lastname'].decode('latin-1')
        i['login'] = i['login'].decode('utf-8')
        i['passwd'] = i['passwd'].decode('utf-8')
        append(i)

    table = new_db['public_interface_members']
    table.insert_many(fixed_items)


def migrate_primers_table(old_db, new_db):
    fixed_items = []
    append = fixed_items.append
    res = old_db.query("select * from primers")
    for i in res:
        del i['id']
        i['gene_code'] = i['geneCode']
        del i['geneCode']
        del i['timestamp']
        if i['primer1'] is None:
            i['primer1'] = ''
        if i['primer2'] is None:
            i['primer2'] = ''
        if i['primer3'] is None:
            i['primer3'] = ''
        if i['primer4'] is None:
            i['primer4'] = ''
        if i['primer5'] is None:
            i['primer5'] = ''
        if i['primer6'] is None:
            i['primer6'] = ''
        append(i)

    table = new_db['public_interface_primers']
    table.insert_many(fixed_items)


def migrate_sequences_table(old_db, new_db):
    fixed_items = []
    append = fixed_items.append
    res = old_db.query("select * from sequences")
    for i in res:
        del i['id']
        i['gene_code'] = i['geneCode']
        del i['geneCode']
        del i['timestamp']

        i['time_created'] = i['dateCreation']
        del i['dateCreation']
        i['time_edited'] = i['dateModification']
        del i['dateModification']
        try:
            i['labPerson'] = i['labPerson'].decode('utf-8')
        except UnicodeDecodeError:
            i['labPerson'] = i['labPerson'].decode('latin-1')
        except AttributeError:
            pass
        append(i)
    table = new_db['public_interface_sequences']
    table.insert_many(fixed_items)


def migrate_taxonsets_table(old_db, new_db):
    fixed_items = []
    append = fixed_items.append
    res = old_db.query("select * from taxonsets")
    for i in res:
        try:
            del i['id']
        except:
            pass
        del i['taxonset_id']

        try:
            i['taxonset_creator'] = i['taxonset_creator'].decode('utf-8')
        except UnicodeDecodeError:
            i['taxonset_creator'] = i['taxonset_creator'].decode('latin-1')
        except AttributeError:
            pass
        append(i)

    table = new_db['public_interface_taxonsets']
    table.insert_many(fixed_items)


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

    def parse_table_genes(self, xml_string):
        this_table = self.tables_prefix + "genes"

        root = ET.fromstring(xml_string)
        for i in root.iter('table_data'):
            if i.attrib['name'] == this_table:
                our_data = i
                break

        for row in our_data.findall('row'):
            id = row.find("./field/[@name='geneCode']")
            print(id.text)





dump_file = sys.argv[1].strip()
with codecs.open(dump_file, "r") as handle:
    dump = handle.read()

tables_prefix = 'voseq_'
parser = ParseXML(dump, tables_prefix)
print(parser.table_genes)
