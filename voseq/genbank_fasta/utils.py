import json

from public_interface.models import Vouchers
from public_interface.models import Sequences
from public_interface.models import Genes


def get_voucher_codes(cleaned_data):
    voucher_codes = []
    if cleaned_data['taxonset'] is not None:
        voucher_codes = json.loads(cleaned_data['taxonset'].taxonset_list)
    if cleaned_data['voucher_codes'] != '':
        voucher_codes += cleaned_data['voucher_codes'].splitlines()
    return set(voucher_codes)


def get_gene_codes(cleaned_data):
    gene_codes = []
    if cleaned_data['geneset'] is not None:
        gene_codes = json.loads(cleaned_data['geneset'].geneset_list)
    if len(cleaned_data['gene_codes']) > 0:
        gene_codes += [i.gene_code for i in cleaned_data['gene_codes']]
    return set(gene_codes)


def get_results(voucher_codes, gene_codes):
    """Returns:
    * List of items with accession numbers that will not be included in
      FASTA or protein files (code, gene_code, accession).
    * FASTA dataset as string.
    * Protein dataset as string.
    """
    items_with_accession = []
    fasta = ''
    for code in voucher_codes:
        v = None
        try:
            v = Vouchers.objects.get(code=code)
        except Vouchers.DoesNotExist:
            continue

        for gene_code in gene_codes:
            s = None
            if v:
                try:
                    s = Sequences.objects.get(code=v, gene_code=gene_code)
                except Sequences.DoesNotExist:
                    continue

                try:
                    g = Genes.objects.get(gene_code=gene_code)
                except Genes.DoesNotExist:
                    continue

                if s.accession.strip() != '':
                    items_with_accession.append((code, gene_code, s.accession))
                else:
                    fasta += '>' + v.genus + '_' + v.species + '_' + code
                    fasta += ' [org=' + v.genus + ' ' + v.species + ']'
                    fasta += ' [Specimen-voucher=' + v.code + ']'
                    fasta += ' [note=' + g.description + ']'
                    fasta += ' [Lineage=]\n'
                    fasta += s.sequences + '\n'
