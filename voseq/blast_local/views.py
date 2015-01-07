from django.shortcuts import render

from .utils import BLAST


def index(request, voucher_code, gene_code):
    blast = BLAST('local', voucher_code, gene_code, test=True)
    blast.save_seqs_to_file()
    blast.create_blast_db()
    blast.save_query_to_file()
    blast.do_blast()
    result = blast.parse_blast_output()
    return render(request, 'blast_local/index.html',
                  {'result': result})
