import re

from django.utils import timezone

from public_interface.models import TaxonSets, GeneSets, Genes
from voseq.celery import app
from create_dataset.models import Dataset
from .utils import CreateDataset


@app.task
def create_dataset(
    taxonset_id, geneset_id, gene_codes_ids, voucher_codes, file_format,
    outgroup, positions, partition_by_positions, translations, aminoacids,
    degen_translations, special, taxon_names, number_genes, introns,
    dataset_obj_id
):
    cleaned_data = dict()
    taxon_set = TaxonSets.objects.get(id=taxonset_id)
    gene_set = GeneSets.objects.get(id=geneset_id)
    cleaned_data['taxonset'] = taxon_set
    cleaned_data['geneset'] = gene_set
    cleaned_data['gene_codes'] = Genes.objects.filter(id__in=gene_codes_ids)
    cleaned_data['voucher_codes'] = voucher_codes
    cleaned_data['file_format'] = file_format
    cleaned_data['outgroup'] = outgroup
    cleaned_data['positions'] = positions
    cleaned_data['partition_by_positions'] = partition_by_positions
    cleaned_data['translations'] = translations
    cleaned_data['aminoacids'] = aminoacids
    cleaned_data['degen_translations'] = degen_translations
    cleaned_data['special'] = special
    cleaned_data['taxon_names'] = taxon_names
    cleaned_data['number_genes'] = number_genes
    cleaned_data['introns'] = introns
    dataset_creator = CreateDataset(cleaned_data)
    dataset_obj = Dataset.objects.get(id=dataset_obj_id)
    dataset_obj.content = dataset_creator.dataset_str
    dataset_obj.completed = timezone.now()
    dataset_obj.save()
    dataset = "{}{}{}".format(
        dataset_creator.dataset_str[0:1500],
        '\n...\n\n\n',
        '#######\nComplete dataset file available for download.\n#######',
    )
    errors = dataset_creator.errors
    warnings = set(dataset_creator.warnings)

    dataset_file_abs = dataset_creator.dataset_file
    if dataset_file_abs is not None:
        dataset_file = re.search(
            '([A-Z]+_[a-z0-9]+\.txt)',
            dataset_file_abs
        ).groups()[0]
    else:
        dataset_file = False

    context = dict()
    context['dataset_file'] = dataset_file
    context['charset_block'] = dataset_creator.charset_block
    context['dataset'] = dataset
    context['dataset_format'] = file_format
    context['errors'] = errors
    context['warnings'] = warnings

