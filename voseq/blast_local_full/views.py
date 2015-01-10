from django.shortcuts import render
from django.conf import settings

from .utils import BLASTFull
from stats.models import Stats


VERSION = settings.VERSION
try:
    STATS = Stats.objects.get(pk=1)
except Stats.DoesNotExist:
    STATS = ''


def index(request, voucher_code, gene_code):
    blast = BLASTFull('full', voucher_code, gene_code)
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
                  }
    )
