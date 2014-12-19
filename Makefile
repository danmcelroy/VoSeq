serve:
	cd voseq; python manage.py runserver --settings=voseq.settings.local

migrations:
	cd voseq; python manage.py makemigrations --settings=voseq.settings.local
	cd voseq; python manage.py migrate --settings=voseq.settings.local

index:
	cd voseq; python manage.py rebuild_index --settings=voseq.settings.local

coverage:
	rm -rf htmlcov .coverage
	coverage run --source voseq voseq/manage.py test -v 2 public_interface --settings=voseq.settings.base
	coverage report -m
	coverage html

test:
	coverage run --source voseq voseq/manage.py test -v 2 public_interface --settings=voseq.settings.base
