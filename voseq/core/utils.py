import re

from django.conf import settings

from stats.models import Stats
from . import exceptions


def get_voucher_codes(cleaned_data):
    """Processes list of voucher codes entered by users.

    It receives data from a **Form class** (`cleaned_data`) and makes sure that
    there are not duplicated voucher codes passed to the dataset builder.
    It also drops voucher codes if specified by users using the double dash:
    `--CP100-10`.

    Args:
        * `form.cleaned_data`: taken from a Form class

    Returns:
        tuple of voucher codes, no dupes, dropped unwanted.
    """
    voucher_codes = tuple()
    if cleaned_data['taxonset'] is not None:
        voucher_codes += tuple(cleaned_data['taxonset'].taxonset_list.splitlines())
    if cleaned_data['voucher_codes'] != '':
        voucher_codes += tuple(cleaned_data['voucher_codes'].splitlines())

    voucher_codes_clean = tuple()
    for i in voucher_codes:
        if re.search('^--', i):
            i_clean = re.sub('^--', '', i)
            voucher_codes_clean += (i_clean,)
        else:
            voucher_codes_clean += (i,)

    voucher_codes_set = tuple()
    for i in voucher_codes_clean:
        if i not in voucher_codes_set and i.strip() != '':
            voucher_codes_set += (i,)

    vouchers_to_drop = []
    for i in voucher_codes:
        if re.search('^--', i):
            vouchers_to_drop.append(re.sub('^--', '', i))

    voucher_codes_filtered = tuple()
    for i in voucher_codes_set:
        if i not in vouchers_to_drop:
            voucher_codes_filtered += (i,)
    return voucher_codes_filtered


def get_gene_codes(cleaned_data):
    """Processes list of gene codes entered by users.

    It receives data from a **Form class** (`cleaned_data`) and makes sure that
    there are not duplicated gene codes passed to the dataset builder.

    Args:
        * `form.cleaned_data`: taken from a Form class

    Returns:
        set of gene codes, sorted, no duplicates.
    """
    gene_codes = []
    if cleaned_data['geneset'] is not None:
        geneset_list = cleaned_data['geneset'].geneset_list.splitlines()
        for i in geneset_list:
            if i not in gene_codes:
                gene_codes.append(i)

    if cleaned_data['gene_codes']:
        for i in cleaned_data['gene_codes']:
            if i.gene_code not in gene_codes:
                gene_codes.append(i.gene_code)

    return tuple(sorted(gene_codes, key=str.lower))


def get_version_stats():
    """Returns version and database statistics for page footer.

    """
    version = settings.VERSION
    try:
        stats = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        stats = ''

    return version, stats


def get_username(request):
    username = 'Guest'
    if request.user.is_authenticated():
        username = request.user.username
    return username


def clean_positions(a_list):
    if 'ALL' in a_list:
        return ['ALL']
    elif '1st' in a_list and '2nd' in a_list and '3rd' in a_list:
        return ['ALL']
    elif '1st' in a_list and '2nd' in a_list and '3rd' not in a_list:
        return ['1st-2nd']
    elif len(a_list) == 2:  # 1st and 3rd, or 2nd and 3rd
        raise exceptions.InadequateCodonPositions(
            "Cannot create dataset for only codon positions {0} and {1}.".format(a_list[0], a_list[1])
        )
    else:
        return a_list
