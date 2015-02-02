from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from public_interface.models import Sequences
from public_interface.models import Vouchers


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
        self.cleaned_data = cleaned_data
        self.dataset_str = self.create_dataset()
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.taxon_names = cleaned_data['taxon_names']

    def create_dataset(self):
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()
        return self.from_seq_objs_to_fasta()

    def create_seq_objs(self):
        """Generate a list of sequence objects. Also takes into account the
        genes passed as geneset.

        Args:
            * ``voucher_codes``: list of vouchers codes, cleaned by our Form.
            * ``gene_codes``: list of gene codes, cleaned by our Form.

        Returns:
            list of sequence objects as produced by BioPython.

        """
        all_seqs = Sequences.objects.all().values('code_id', 'gene_code', 'sequences').order_by('code_id')
        for s in all_seqs:
            code = s['code_id'].lower()
            gene_code = s['gene_code'].lower()
            if code in self.voucher_codes and gene_code in self.gene_codes:
                seq = Seq(s['sequences'])
                seq_obj = SeqRecord(seq)
                seq_obj.id = code
                seq_obj.name = gene_code

                if gene_code not in self.seq_objs:
                    self.seq_objs[gene_code] = []
                self.seq_objs[gene_code].append(seq_obj)

    def from_seq_objs_to_fasta(self):
        """Take a list of BioPython's sequence objects and return a FASTA string

        Returns:
            The FASTA string will contain the gene_code at the top as if it were
            another FASTA gene sequence.

        """
        fasta_str = []
        append = fasta_str.append

        for gene_code in self.seq_objs:
            this_gene = None
            for seq_record in self.seq_objs[gene_code]:
                if this_gene is None:
                    this_gene = seq_record.name
                    seq_str = '>' + this_gene + '\n' + '--------------------'
                    append(seq_str)
                if this_gene != seq_record.name:
                    this_gene = seq_record.name
                    seq_str = '>' + this_gene + '\n' + '--------------------'
                    append(seq_str)
                seq_str = '>' + seq_record.id.upper() + '\n' + str(seq_record.seq)
                append(seq_str)

        return '\n'.join(fasta_str)

    def get_taxon_names_for_taxa(self):
        """Returns list of dicts: {'taxon': 'name'}

        Takes list of voucher_codes and list of taxon_names from cleaned form.

        Args:
            * ``voucher_codes``
            * ``taxon_names``

        Returns:
            List of dictionaries with data.

        """
        vouchers_with_taxon_names = []
        append = vouchers_with_taxon_names.append

        all_vouchers = Vouchers.objects.all().order_by('code').values()
        for voucher in all_vouchers:
            code = voucher['code'].lower()
            if code in self.voucher_codes:
                obj = dict()
                for taxon_name in self.taxon_names:
                    taxon_name = taxon_name.lower()
                    obj[taxon_name] = voucher[taxon_name]
                append(obj)

        return vouchers_with_taxon_names
