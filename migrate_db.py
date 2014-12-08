#!-*- encoding: utf-8 -*-
"""
Needs an XML file with the database dump from MySQL:

> mysqldump --xml database > dump.xml
"""
import datetime
import dataset
import json

with open("config.json", "r") as f:
    settings = json.loads(f.read())


new_db_url = 'postgresql://' + settings['POSTGRES_DB_USER'] + ':' + settings['POSTGRES_DB_PASS'] + '@' \
             + 'localhost:5432/' + settings['POSTGRES_DB_NAME']
old_db = dataset.connect(old_db_url)
new_db = dataset.connect(new_db_url)


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



migrate_sequences_table(old_db, new_db)
