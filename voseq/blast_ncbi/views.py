from django.shortcuts import render

from .utils import BLASTNcbi


def index(request, voucher_code, gene_code):
    blast = BLASTNcbi(voucher_code, gene_code)
    blast.save_query_to_file()

    blast.do_blast()
    result = blast.parse_blast_output()
    blast.delete_query_output_files()
    return render(request, 'blast_local/index.html',
                  {'result': result})
