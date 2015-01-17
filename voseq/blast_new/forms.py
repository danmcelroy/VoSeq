from django import forms

from public_interface.models import Genes


class BLASTNewForm(forms.Form):
    name = forms.CharField(
        label='Name',
        max_length=100,
    )
    sequence = forms.CharField(
        label='Sequence',
        required=True,
        widget=forms.Textarea,
    )

    queryset = Genes.objects.all()
    CHOICES = [(i.gene_code, i.gene_code) for i in queryset]
    gene_codes = forms.ChoiceField(
        choices=CHOICES,
        label='Genes',
        widget=forms.CheckboxSelectMultiple(),
    )
