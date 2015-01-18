from django import forms

from public_interface.models import Genes
from public_interface.models import GeneSets
from public_interface.models import TaxonSets


class GenBankFastaForm(forms.Form):
    taxonset = forms.ModelChoiceField(
        TaxonSets.objects.all(),
        label='Choose taxonset',
        required=False,
        empty_label='Choose taxonset',
    )
    geneset = forms.ModelChoiceField(
        GeneSets.objects.all(),
        label='Choose geneset',
        required=False,
        empty_label='Choose geneset',
    )
    gene_codes = forms.ModelMultipleChoiceField(
        Genes.objects.all(),
        label='Check yo select your alignment/gene',
        required=False,
        widget=forms.CheckboxSelectMultiple(),
    )
    voucher_codes = forms.CharField(
        widget=forms.Textarea,
        label='... and/or a list of voucher codes',
        required=False,
    )

    def clean(self):
        """Overwriting validator method of class form."""
        cleaned_data = super(GenBankFastaForm, self).clean()
        print(cleaned_data)
        taxonset = cleaned_data.get("taxonset")
        voucher_codes = cleaned_data.get("voucher_codes")

        if taxonset is None and voucher_codes.strip() == '':
            raise forms.ValidationError("You need to enter at least some "
                                        "voucher codes or select a taxonset.")
