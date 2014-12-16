serve:
	cd voseq; python manage.py runserver --settings=voseq.settings.local

migrations:
	cd voseq; python manage.py makemigrations --settings=voseq.settings.local
	cd voseq; python manage.py migrate --settings=voseq.settings.local

