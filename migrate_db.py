# This needs python2.7
import datetime
import dataset
import json

with open("config.json", "r") as f:
    settings = json.loads(f.read())


old_db_url = 'mysql://' + settings['DB_USER'] + ':' + settings['DB_PASS'] + '@' + settings['DB_HOST'] + '/' + settings['DB_NAME']
new_db_url = 'postgresql://' + settings['POSTGRES_DB_USER'] + ':' + settings['POSTGRES_DB_PASS'] + '@' \
             + 'localhost:5432/' + settings['POSTGRES_DB_NAME']
old_db = dataset.connect(old_db_url)
new_db = dataset.connect(new_db_url)

fixed_items = []
append = fixed_items.append

print(old_db['genes'].columns)
res = old_db.query("select * from genes")
for i in res:
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
