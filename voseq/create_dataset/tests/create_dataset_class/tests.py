from django.test import TestCase

from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import GeneSets
from public_interface.models import Vouchers
from public_interface.models import TaxonSets
from create_dataset.utils import CreateDataset


class TestCreateDataset(TestCase):
    def setUp(self):
        """
        genes = []
        for i in ['abc1', 'BC1', 'ABC2', 'CC', 'xaz', 'XYZ']:
            genes.append(Genes(gene_code=i))
        Genes.objects.bulk_create(genes)
        """
        genes = ['abc1', 'BC1', 'ABC2', 'CC', 'xaz', 'XYZ']
        GeneSets(geneset_name='6genes', geneset_creator='Carlos Pena',
                 geneset_list="\n".join(genes)).save()

        taxa = ['CP100-10', 'CP100-11', 'CP100-12']
        TaxonSets(taxonset_name='all_taxa', taxonset_creator='Carlos Pena',
                  taxonset_list="\n".join(taxa)).save()

        tmp = []
        for taxon in taxa:
            tmp.append(Vouchers(code=taxon))
        Vouchers.objects.bulk_create(tmp)

        tmp = []
        for gene in genes:
            tmp.append(Genes(gene_code=gene, genetic_code=1, length=600,
                             reading_frame=1))
        Genes.objects.bulk_create(tmp)

        for gene in genes:
            for taxon in taxa:
                voucher = Vouchers.objects.get(code=taxon)
                Sequences(code=voucher, gene_code=gene, sequences='ACTG',
                          genbank=False).save()

        self.cleaned_data = {
            'voucher_codes': '',
            'positions': ['ALL'],
            'introns': 'YES',
            'gene_codes': [],
            'partition_by_positions': 'ONE',
            'number_genes': None,
            'geneset': GeneSets.objects.get(geneset_name='6genes'),
            'file_format': 'PHY',
            'taxon_names': ['CODE', 'GENUS', 'SPECIES'],
            'special': False,
            'taxonset': TaxonSets.objects.get(taxonset_name='all_taxa'),
            'outgroup': '',
            'aminoacids': False,
            'translations': False,
            'degen_translations': 'NORMAL'
        }

    def test_seq_objs_have_sorted_gene_codes(self):
        result = CreateDataset(self.cleaned_data)
        self.assertEqual(['abc1', 'ABC2', 'BC1', 'CC', 'xaz', 'XYZ'],
                         list(result.seq_objs))
