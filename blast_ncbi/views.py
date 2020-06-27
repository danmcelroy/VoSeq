from django.http import HttpRequest, HttpResponse
from django.shortcuts import render

from core.utils import get_context
from .utils import BLASTNcbi


def index(request: HttpRequest, voucher_code: str, gene_code: str) -> HttpResponse:
    """Executes blast of sequence against NCBI genbank"""
    blast = BLASTNcbi(blast_type="remote", voucher_code=voucher_code,
                      gene_code=gene_code)
    blast.save_query_to_file()

    blast.do_blast()
    result = blast.parse_blast_output()
    blast.delete_query_output_files()
    context = get_context(request)
    context["result"] = result
    return render(request, 'blast_local/index.html', context)
