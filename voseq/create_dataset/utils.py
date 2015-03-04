from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import flatten_taxon_names_dict
from public_interface.models import Sequences
from public_interface.models import Vouchers
from .dataset import Dataset


class CreateFasta(Dataset):
    pass


class CreateTNT(Dataset):
    def __init__(self, *args, **kwargs):
        super(CreateTNT, self).__init__(*args, **kwargs)
        self.number_taxa = None
        self.number_chars = None
        print(self.seq_objs)


class CreateDataset(object):
    """
    Accept form input to create a dataset in several formats, codon positions,
    for list of codes and genes. Also takes into account the vouchers passed as
    taxonset.

    Attributes:
        ``dataset_str``: output dataset to pass to users.

    """
    def __init__(self, cleaned_data):
        print(">>>>>>_init", cleaned_data)
        self.errors = []
        self.seq_objs = dict()
        self.codon_positions = cleaned_data['positions']
        self.file_format = cleaned_data['file_format']
        self.partition_by_positions = cleaned_data['partition_by_positions']
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.taxon_names = cleaned_data['taxon_names']
        self.warnings = []
        self.dataset_str = self.create_dataset()

    def create_dataset(self):
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

        if self.file_format == 'FASTA':
            fasta = CreateFasta(self.codon_positions, self.partition_by_positions, self.seq_objs, self.gene_codes, self.file_format)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            return fasta_dataset

        if self.file_format == 'TNT':
            tnt = CreateTNT(self.codon_positions, self.partition_by_positions, self.seq_objs, self.gene_codes, self.file_format)
            tnt_dataset = tnt.from_seq_objs_to_dataset()
            self.warnings += tnt.warnings
            return tnt_dataset

    def create_seq_objs(self):
        """Generate a list of sequence objects. Also takes into account the
        genes passed as geneset.

        Returns:
            list of sequence objects as produced by BioPython.

        """
        our_taxon_names = self.get_taxon_names_for_taxa()

        all_seqs = Sequences.objects.all().values('code_id', 'gene_code', 'sequences').order_by('code_id')
        for s in all_seqs:
            code = s['code_id'].lower()
            gene_code = s['gene_code'].lower()
            if code in self.voucher_codes and gene_code in self.gene_codes:
                seq = Seq(s['sequences'])
                seq_obj = SeqRecord(seq)
                seq_obj.id = flatten_taxon_names_dict(our_taxon_names[code])
                if 'GENECODE' in self.taxon_names:
                    seq_obj.id += '_' + gene_code
                seq_obj.name = gene_code

                if gene_code not in self.seq_objs:
                    self.seq_objs[gene_code] = []
                self.seq_objs[gene_code].append(seq_obj)

    def get_taxon_names_for_taxa(self):
        """Returns dict: {'CP100-10': {'taxon': 'name'}}

        Takes list of voucher_codes and list of taxon_names from cleaned form.

        Returns:
            Dictionary with data, also as dicts.

        """
        vouchers_with_taxon_names = {}

        all_vouchers = Vouchers.objects.all().order_by('code').values('code', 'orden', 'superfamily',
                                                                      'family', 'subfamily', 'tribe',
                                                                      'subtribe', 'genus', 'species',
                                                                      'subspecies', 'author', 'hostorg',)
        for voucher in all_vouchers:
            code = voucher['code'].lower()
            if code in self.voucher_codes:
                obj = dict()
                for taxon_name in self.taxon_names:
                    if taxon_name != 'GENECODE':
                        taxon_name = taxon_name.lower()
                        obj[taxon_name] = voucher[taxon_name]
                vouchers_with_taxon_names[code] = obj

        return vouchers_with_taxon_names
