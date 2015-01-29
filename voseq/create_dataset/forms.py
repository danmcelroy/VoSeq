from django import forms

from core.forms import BaseDatasetForm


class CreateDatasetForm(BaseDatasetForm):
    file_format = forms.ChoiceField(
        label='Choose file format',
        choices=[
            ('TNT', 'TNT format'),
            ('NEXUS', 'NEXUS format'),
            ('PHY', 'PHYLIP format'),
            ('FASTA', 'Unaligned FASTA format'),
        ],
        widget=forms.RadioSelect(),
        required=True,
        initial='FASTA',
    )

    outgroup = forms.CharField(
        label='Outgroup (code, for NEXUS and TNT)',
        help_text='Voucher code for using that specimen\'s sequence as '
                  'outgroup in NEXUS and TNT datasets.',
        required=False,
    )

    positions = forms.ChoiceField(
        label='Positions',
        help_text='codon positions',
        choices=[
            ('ALL', 'all'),
            ('1st', '1st'),
            ('2nd', '2nd'),
            ('3rd', '3rd'),
        ],
        widget=forms.RadioSelect(),
        initial='ALL',
        required=True,
    )

    partition_by_positions = forms.ChoiceField(
        label='Partition by positions',
        choices=[
            ('ONE', 'as one'),
            ('EACH', 'each'),
            ('1st-3rd', '1st-2nd, 3rd'),
        ],
        widget=forms.RadioSelect(),
        initial='ONE',
        required=True,
    )

    translations = forms.BooleanField(
        label='Degen(erated)',
        widget=forms.CheckboxInput(),
        required=False,
    )

    aminoacids = forms.BooleanField(
        label='Amino acids',
        widget=forms.CheckboxInput(),
        required=False,
    )

    degen_translations = forms.ChoiceField(
        label='Degenerated translations)',
        choices=[
            ('NORMAL', 'normal'),
            ('S', 'S'),
            ('Z', 'Z'),
            ('SZ', 'SZ'),
        ],
        widget=forms.RadioSelect(),
        initial='NORMAL',
        required=True,
    )

    special = forms.BooleanField(
        label='Special',
        widget=forms.CheckboxInput(),
        required=False,
    )

    taxon_names = forms.MultipleChoiceField(
        label='What info do you want in the taxon names?',
        choices=[
            ('CODE', 'Code'),
            ('ORDER', 'Order'),
            ('FAMILY', 'Family'),
            ('SUBFAMILY', 'Subfamily'),
            ('TRIBE', 'Tribe'),
            ('SUBTRIBE', 'Subtribe'),
            ('GENUS', 'Genus'),
            ('SPECIES', 'Species'),
            ('SUBSPECIES', 'Subspecies'),
            ('AUCTOR', 'Auctor'),
            ('HOSTORG', 'Host org.'),
            ('GENECODE', 'Gene code'),
        ],
        widget=forms.CheckboxSelectMultiple(),
        initial=['CODE', 'GENUS', 'SPECIES'],
        required=False,
        help_text='If taxon_names is None, use standart code_genus_species',
    )

    exclude = forms.ChoiceField(
        help_text='If dataset is for single gene, exclude taxa missing this gene? '
                  'Otherwise include taxon with ? as sequences.',
        choices=[
            ('YES', 'yes'),
            ('NO', 'no'),
        ],
        widget=forms.RadioSelect(),
        initial='YES',
        required=True,
    )

    number_genes = forms.IntegerField(
        label='Minimum number of genes',
        required=False,
    )

    introns = forms.ChoiceField(
        choices=[
            ('YES', 'yes'),
            ('NO', 'no'),
        ],
        widget=forms.RadioSelect(),
        initial='YES',
        required=True,
        help_text='Ignore introns?',
    )
