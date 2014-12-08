serve:
	cd genomeseq; python manage.py runserver --settings=genomeseq.settings.local

migrations:
	cd genomeseq; python manage.py makemigrations --settings=genomeseq.settings.local
	cd genomeseq; python manage.py migrate --settings=genomeseq.settings.local
