from django.contrib.auth.decorators import login_required
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render

from core.utils import get_context
from .utils import BLASTFull


@login_required
def index(request: HttpRequest, voucher_code: str, gene_code: str
          ) -> HttpResponse:
    """Execute a blast of sequence against all sequences in the database"""
    blast = BLASTFull('full', voucher_code, gene_code)
    blast.save_seqs_to_file()

    if not blast.is_blast_db_up_to_date():
        blast.create_blast_db()

    blast.save_query_to_file()
    blast.do_blast()
    result = blast.parse_blast_output()
    blast.delete_query_output_files()
    context = get_context(request)
    context['result'] = result
    return render(request, 'blast_local/index.html', context)
