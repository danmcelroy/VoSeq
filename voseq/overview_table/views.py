from django.shortcuts import render

from core.utils import get_version_stats
from core.utils import get_username
from public_interface.models import Sequences
from public_interface.models import Vouchers


def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    vouchers = Vouchers.objects.all().values('code', 'orden', 'family',
                                             'subfamily', 'genus', 'species',
                                             'hostorg').order_by('code')
    vouchers = convert_to_dict(vouchers)
    sequences = Sequences.objects.all().values('code', 'gene_code', 'sequences').order_by('code')

    for sequence in sequences:
        code = sequence['code']
        trimmed_sequence = {
            'length': len(sequence['sequences']),
            'gene_code': sequence['gene_code'],
        }
        this_dict = vouchers[code]
        this_dict['sequences'].append(trimmed_sequence)
        vouchers[code] = this_dict

    return render(request,
                  'overview_table/index.html',
                  {
                      'username': username,
                      'version': version,
                      'stats': stats,
                      'data': vouchers,
                  },
                  )


def convert_to_dict(vouchers):
    out = dict()
    for voucher in vouchers:
        voucher['sequences'] = []
        code = voucher['code']
        out[code] = voucher
    return out


def join_dictionaries(x, y):
    return dict(list(x.items()) + list(y.items()))
