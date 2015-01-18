from django import forms

from public_interface.models import Genes
from public_interface.models import GeneSets
from public_interface.models import TaxonSets


class GenBankFastaForm(forms.Form):
    taxonset = forms.ModelChoiceField(
        TaxonSets.objects.all(),
        label='Choose taxonset',
        required=False,
    )
    geneset = forms.ModelChoiceField(
        GeneSets.objects.all(),
        label='Choose geneset',
        required=False,
    )
    gene_codes = forms.ModelMultipleChoiceField(
        Genes.objects.all(),
        label='Check yo select your alignment/gene',
        required=False,
        widget=forms.CheckboxSelectMultiple(),
    )
    gene_code_list = forms.CharField(
        widget=forms.Textarea,
        label='... and/or a list of codes',
        required=False,
    )
