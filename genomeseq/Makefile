
coverage:
	rm -rf htmlcov
	coverage run --source=. manage.py test -v 2 public_interface --settings=genomeseq.settings.base
	coverage report -m
	coverage html