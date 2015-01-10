from django.conf import settings
from django.shortcuts import render

from .utils import BLAST
from stats.models import Stats


VERSION = settings.VERSION
STATS = Stats.objects.get(pk=1)
print("AAA", STATS)


def index(request, voucher_code, gene_code):
    blast = BLAST('local', voucher_code, gene_code)
    blast.save_seqs_to_file()

    if blast.is_blast_db_up_to_date() is False:
        blast.create_blast_db()

    blast.save_query_to_file()
    blast.do_blast()
    result = blast.parse_blast_output()
    blast.delete_query_output_files()
    return render(request, 'blast_local/index.html',
                  {
                      'result': result,
                      'version': VERSION,
                      'stats': STATS,
                  },
                  )
