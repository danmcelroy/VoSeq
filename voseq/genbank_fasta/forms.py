from django import forms

from public_interface.models import Gene
from public_interface.models import GeneSet
from public_interface.models import TaxonSet


class GenBankFastaForm(forms.Form):
    taxonset = forms.ModelChoiceField(
        TaxonSet.objects.all(),
        label='Choose taxonset',
        required=False,
        empty_label='Choose taxonset',
    )
    geneset = forms.ModelChoiceField(
        GeneSet.objects.all(),
        label='Choose geneset',
        required=False,
        empty_label='Choose geneset',
    )
    gene_codes = forms.ModelMultipleChoiceField(
        Gene.objects.all(),
        label='Check to select your alignment/gene',
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
        taxonset = cleaned_data.get("taxonset")
        voucher_codes = cleaned_data.get("voucher_codes")

        geneset = cleaned_data.get("geneset")
        gene_codes = cleaned_data.get("gene_codes")
        print(gene_codes)

        if taxonset is None and voucher_codes.strip() == '':
            raise forms.ValidationError("You need to enter at least some "
                                        "voucher codes or select a taxonset.")

        if geneset is None and len(gene_codes) < 1:
            raise forms.ValidationError("You need to enter at least some "
                                        "gene codes or select a geneset.")
