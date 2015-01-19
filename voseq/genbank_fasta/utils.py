import json


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
