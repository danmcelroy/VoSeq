from django import forms

from core.forms import BaseDatasetForm


class VoucherTableForm(BaseDatasetForm):
    voucher_info = forms.MultipleChoiceField(
        label='Voucher info',
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
        label='Locality and Collector info',
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
