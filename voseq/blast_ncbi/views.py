from django.shortcuts import render
from django.conf import settings

from .utils import BLASTNcbi
from stats.models import Stats


def index(request, voucher_code, gene_code):
    VERSION = settings.VERSION
    try:
        STATS = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        STATS = ''

    blast = BLASTNcbi(voucher_code, gene_code)
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
