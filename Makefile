serve: index stats
	cd voseq; python manage.py runserver --settings=voseq.settings.local

migrations:
	cd voseq; python manage.py makemigrations --settings=voseq.settings.local
	cd voseq; python manage.py migrate --settings=voseq.settings.local

import:
	python voseq/manage.py migrate_db --dumpfile=test_db_dump.xml --settings=voseq.settings.local

index:
	cd voseq; python manage.py rebuild_index --settings=voseq.settings.local

stats:
	cd voseq; python manage.py create_stats --settings=voseq.settings.local

coverage:
	rm -rf htmlcov .coverage
	coverage run --source voseq voseq/manage.py test -v 2 public_interface blast_local blast_local_full blast_ncbi --settings=voseq.settings.base
	coverage report -m
	coverage html

test:
	coverage run --source voseq voseq/manage.py test -v 2 public_interface blast_local blast_local_full blast_ncbi --settings=voseq.settings.base
