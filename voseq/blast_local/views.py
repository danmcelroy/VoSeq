import logging
from subprocess import CalledProcessError

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render

from core.utils import get_version_stats, get_username, BLAST


log = logging.getLogger(__name__)


def index(request: HttpRequest, voucher_code: str, gene_code: str) -> HttpResponse:
    """Runs a local blast for voucher code and gene code

    Show results to user in a table
    """
    version, stats = get_version_stats()
    username = get_username(request)

    blast = BLAST('local', voucher_code, gene_code)
    blast.save_seqs_to_file()

    if not blast.is_blast_db_up_to_date():
        try:
            blast.create_blast_db()
        except CalledProcessError:
            log.warning("there are no sequences for gene %s", gene_code)

    was_sequence_saved = blast.save_query_to_file()
    if was_sequence_saved:
        blast.do_blast()
        result = blast.parse_blast_output()
        blast.delete_query_output_files()
    else:
        result = {
            "error": "Query sequence has no valid codons, only question marks",
        }
    return render(
        request, 'blast_local/index.html',
        {
            'username': username,
            'result': result,
            'version': version,
            'stats': stats,
        },
    )
