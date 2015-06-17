from collections import OrderedDict

from django.shortcuts import render

from core.utils import get_version_stats
from core.utils import get_username
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers


def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    vouchers_with_sequences = get_vouchers_with_seqs()
    vouchers_without_sequences = get_vouchers_without_seqs()
    # vouchers = vouchers_without_sequences + vouchers_with_sequences

    genes = Genes.objects.all().order_by('gene_code')

    return render(request,
                  'overview_table/index.html',
                  {
                      'username': username,
                      'version': version,
                      'stats': stats,
                      'data': vouchers_with_sequences,
                      'genes': genes,
                  },
                  )


def get_vouchers_with_seqs():
    seqs = Sequences.objects.select_related('code').order_by('code').values('code', 'gene_code')
    return seqs


def get_vouchers_without_seqs():
    vouchers = Vouchers.objects.all().filter(sequences__code__isnull=True).values('code', 'orden', 'family', 'subfamily', 'genus', 'species', 'hostorg').order_by('code')
    return vouchers


def get_empty_sequences():
    genes = Genes.objects.all().order_by('gene_code')
    empty_sequences = OrderedDict()
    for gene in genes:
        empty_sequences[gene.gene_code] = ''
    return empty_sequences


def join_dictionaries(x, y):
    return dict(list(x.items()) + list(y.items()))
