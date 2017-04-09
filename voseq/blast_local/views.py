from django.shortcuts import render

from core.utils import get_version_stats
from core.utils import get_username
from core.utils import BLAST


def index(request, voucher_code, gene_code):
    version, stats = get_version_stats()
    username = get_username(request)

    blast = BLAST('local', voucher_code, gene_code)
    blast.save_seqs_to_file()

    if not blast.is_blast_db_up_to_date():
        blast.create_blast_db()

    good_sequence = blast.save_query_to_file()
    if good_sequence:
        blast.do_blast()
        result = blast.parse_blast_output()
        blast.delete_query_output_files()
    else:
        result = {"error": "Query sequence has no valid codons, only question marks"}
    return render(
        request, 'blast_local/index.html',
        {
            'username': username,
            'result': result,
            'version': version,
            'stats': stats,
        },
    )
