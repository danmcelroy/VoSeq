from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import flatten_taxon_names_dict
from .dataset import CreateGenbankFasta
from .dataset import CreateFasta
from .dataset import CreatePhylip
from .dataset import CreateNEXUS
from .dataset import CreateTNT
from public_interface.models import Genes
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
        self.errors = []
        self.seq_objs = dict()
        self.minimum_number_of_genes = cleaned_data['number_genes']
        self.aminoacids = cleaned_data['aminoacids']
        self.codon_positions = cleaned_data['positions']
        self.file_format = cleaned_data['file_format']
        self.partition_by_positions = cleaned_data['partition_by_positions']
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.gene_codes_and_lengths = None
        self.taxon_names = cleaned_data['taxon_names']
        self.voucher_codes_metadata = dict()
        self.gene_codes_metadata = self.get_gene_codes_metadata()
        self.warnings = []
        self.outgroup = cleaned_data['outgroup']
        self.dataset_file = None
        self.aa_dataset_file = None
        self.charset_block = None
        self.dataset_str = self.create_dataset()

    def create_dataset(self):
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

        if self.file_format == 'GenbankFASTA':
            fasta = CreateGenbankFasta(self.codon_positions, self.partition_by_positions,
                                       self.seq_objs, self.gene_codes, self.voucher_codes,
                                       self.file_format)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            self.dataset_file = fasta.dataset_file
            self.aa_dataset_file = fasta.aa_dataset_file
            return fasta_dataset

        if self.file_format == 'FASTA':
            fasta = CreateFasta(self.codon_positions, self.partition_by_positions,
                                self.seq_objs, self.gene_codes, self.voucher_codes,
                                self.file_format)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            self.dataset_file = fasta.dataset_file
            return fasta_dataset

        if self.file_format == 'PHY':
            phy = CreatePhylip(self.codon_positions, self.partition_by_positions,
                               self.seq_objs, self.gene_codes, self.voucher_codes,
                               self.file_format, self.outgroup, self.voucher_codes_metadata,
                               self.minimum_number_of_genes, self.aminoacids)
            phylip_dataset = phy.from_seq_objs_to_dataset()
            self.warnings += phy.warnings
            self.dataset_file = phy.dataset_file
            self.charset_block = phy.charset_block
            return phylip_dataset

        if self.file_format == 'TNT':
            tnt = CreateTNT(self.codon_positions, self.partition_by_positions,
                            self.seq_objs, self.gene_codes, self.voucher_codes,
                            self.file_format, self.outgroup, self.voucher_codes_metadata,
                            self.minimum_number_of_genes, self.aminoacids)
            tnt_dataset = tnt.from_seq_objs_to_dataset()
            self.warnings += tnt.warnings
            self.dataset_file = tnt.dataset_file
            return tnt_dataset

        if self.file_format == 'NEXUS':
            nexus = CreateNEXUS(self.codon_positions, self.partition_by_positions,
                                self.seq_objs, self.gene_codes, self.voucher_codes,
                                self.file_format, self.outgroup, self.voucher_codes_metadata,
                                self.minimum_number_of_genes, self.aminoacids)
            nexus_dataset = nexus.from_seq_objs_to_dataset()
            self.warnings += nexus.warnings
            self.dataset_file = nexus.dataset_file
            return nexus_dataset

    def create_seq_objs(self):
        """Generate a list of sequence objects. Also takes into account the
        genes passed as geneset.

        Returns:
            list of sequence objects as produced by BioPython.
        """
        # We might need to update our list of vouches and genes
        gene_codes = set()

        our_taxon_names = self.get_taxon_names_for_taxa()

        all_seqs = self.get_all_sequences()

        for code in self.voucher_codes:
            for gene_code in self.gene_codes:
                if gene_code in all_seqs[code]:
                    seq_obj = self.create_seq_record(all_seqs[code][gene_code])
                    seq_obj.id = flatten_taxon_names_dict(our_taxon_names[code])
                    if 'GENECODE' in self.taxon_names:
                        seq_obj.id += '_' + gene_code
                    seq_obj.name = gene_code
                    seq_obj.description = code
                    self.voucher_codes_metadata[code] = seq_obj.id

                    if gene_code not in self.seq_objs:
                        self.seq_objs[gene_code] = tuple()
                    self.seq_objs[gene_code] += (seq_obj,)
                else:
                    print(">>> there is no seq for gene_code %s and code %s" % (gene_code, code))
                    seq = Seq('?' * self.gene_codes_metadata[gene_code])
                    empty_seq_obj = SeqRecord(seq)
                    try:
                        empty_seq_obj.id = self.voucher_codes_metadata[code]
                    except KeyError:
                        pass
                    empty_seq_obj.name = gene_code
                    empty_seq_obj.description = code

                    if gene_code not in self.seq_objs:
                        self.seq_objs[gene_code] = tuple()
                    self.seq_objs[gene_code] += (empty_seq_obj,)

        self.gene_codes = list(gene_codes)
        self.add_missing_seqs()

    def get_all_sequences(self):
        # Return sequences as dict of lists containing sequence and related data
        seqs_dict = {}

        all_seqs = Sequences.objects.all().values('code_id',
                                                  'gene_code',
                                                  'sequences').order_by('code_id')
        for seq in all_seqs:
            code = seq['code_id']
            gene_code = seq['gene_code']

            if code in self.voucher_codes and gene_code in self.gene_codes:
                if code not in seqs_dict:
                    seqs_dict[code] = {gene_code: ''}
                seqs_dict[code][gene_code] = seq

        """
        seqs_with_old_order = tuple()
        for code in self.voucher_codes:
            try:
                seqs_with_old_order += (seqs_dict[code],)
            except KeyError:
                self.warnings += ['Could not find sequences for voucher %s' % code]
        """
        return seqs_dict

    def create_seq_record(self, s):
        """
        Adds ? if the sequence is not long enough
        :param s:
        :return:
        """
        gene_code = s['gene_code']
        length = self.gene_codes_metadata[gene_code]
        sequence = s['sequences']
        length_difference = length - len(sequence)

        sequence += '?' * length_difference
        seq = Seq(sequence)
        seq_obj = SeqRecord(seq)
        return seq_obj

    def add_missing_seqs(self):
        """
        Loops over the created seq_objects and adds sequences full of ? if
        those where not found in our database.

        Uses the updated lists of voucher_codes and gene_codes
        """
        new_seq_objs = OrderedDict()

        for gene_code in self.seq_objs:
            new_seq_objs[gene_code] = tuple()

            i = 0
            for seq_obj in self.seq_objs[gene_code]:
                code = seq_obj.description
                if code == self.voucher_codes[i]:
                    new_seq_objs[gene_code] += (seq_obj,)
                    i += 1
                else:
                    seq = Seq('?' * self.gene_codes_metadata[gene_code])
                    empty_seq_obj = SeqRecord(seq)
                    try:
                        empty_seq_obj.id = self.voucher_codes_metadata[self.voucher_codes[i]]
                    except KeyError:
                        pass
                    empty_seq_obj.name = gene_code
                    empty_seq_obj.description = self.voucher_codes[i]
                    new_seq_objs[gene_code] += (empty_seq_obj,)
                    i += 1
        self.seq_objs = new_seq_objs

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
            code = voucher['code']
            if code in self.voucher_codes:
                obj = dict()
                for taxon_name in self.taxon_names:
                    if taxon_name != 'GENECODE':
                        taxon_name = taxon_name.lower()
                        obj[taxon_name] = voucher[taxon_name]
                vouchers_with_taxon_names[code] = obj

        return vouchers_with_taxon_names

    def get_gene_codes_metadata(self):
        """
        :return: dictionary with genecode and base pair number.
        """
        queryset = Genes.objects.all().values('gene_code', 'length')
        gene_codes_metadata = dict()
        for i in queryset:
            gene_code = i['gene_code']
            gene_codes_metadata[gene_code] = i['length']
        return gene_codes_metadata
