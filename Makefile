admin:
	python voseq/manage.py createsuperuser

serve:  index stats
	python voseq/manage.py runserver

migrations:
	python voseq/manage.py makemigrations
	python voseq/manage.py migrate

import:
	python voseq/manage.py migrate_db --dumpfile=test_db_dump.xml

index:
	python voseq/manage.py rebuild_index

stats:
	python voseq/manage.py create_stats

coverage: test
	coverage report -m
	coverage html

test:
	rm -rf htmlcov .coverage
	coverage run --source voseq voseq/manage.py test -v 2 core public_interface blast_local blast_local_full blast_ncbi blast_new stats view_genes genbank_fasta