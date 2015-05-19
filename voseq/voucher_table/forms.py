from django import forms

from core.forms import BaseDatasetForm


class VoucherTableForm(BaseDatasetForm):
    voucher_info = forms.MultipleChoiceField(
        label='Voucher info:',
        choices=[
            ('CODE', 'Code'),
            ('ORDEN', 'Order'),
            ('SUPERFAMILY', 'Superfamily'),
            ('FAMILY', 'Family'),
            ('SUBFAMILY', 'Subfamily'),
            ('TRIBE', 'Tribe'),
            ('SUBTRIBE', 'Subtribe'),
            ('GENUS', 'Genus'),
            ('SPECIES', 'Species'),
            ('SUBSPECIES', 'Subspecies'),
            ('AUTHOR', 'Author'),
            ('HOSTORG', 'Host org.'),
        ],
        widget=forms.CheckboxSelectMultiple(),
        initial=['CODE', 'GENUS', 'SPECIES'],
        required=False,
        help_text='If taxon_names is None, use standart code_genus_species',
    )
    collector_info = forms.MultipleChoiceField(
        label='Locality and collector info:',
        choices=[
            ('COUNTRY', 'Country'),
            ('LOCALITY', 'Locality'),
            ('COLLECTOR', 'Collector'),
            ('COLL. DATE', 'Coll. date'),
            ('DETERMINED BY', 'Determined by'),
            ('ALTITUDE', 'Altitude'),
            ('LATITUDE', 'Latitude'),
            ('LONGITUDE', 'Longitude'),
        ],
        widget=forms.CheckboxSelectMultiple(),
        required=False,
    )
    gene_info = forms.ChoiceField(
        label='Choose what gene info to display:',
        choices=[
            ('NUMBER OF BASES', 'Number of bases'),
            ('ACCESSION NUMBER', 'Accession number'),
            ('EXIST OR EMPTY', 'X/- (exists/empty)'),
        ],
        widget=forms.RadioSelect(),
        required=False,
    )
    field_delimitor = forms.ChoiceField(
        label='Choose your field delimitor:',
        choices=[
            ('COMMA', 'comma'),
            ('TAB', 'tab'),
        ],
        widget=forms.RadioSelect(),
        required=False,
    )
