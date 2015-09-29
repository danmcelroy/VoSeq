from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from seqrecord_expanded import SeqRecordExpanded
from dataset_creator import Dataset

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import clean_positions
from core.utils import flatten_taxon_names_dict
from .dataset import CreateGenbankFasta
from .dataset import CreateFasta
from .dataset import CreateTNT
from .mega import CreateMEGA
from .nexus import DatasetHandler
from .phylip import CreatePhylip
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers


class CreateDataset(object):
    """
    Accepts form input to create a dataset in several formats, codon positions,
    for list of codes and genes. Also takes into account the vouchers passed as
    taxonset.

    Attributes:
        ``codon_positions``: For now is a list. It is cleaned to avoid redundant
                      information such as having ['ALL', '1st'] in there as 'ALL'
                      overrides '1st'.
        ``seq_objs``: Ordered_dict by gene_code. Keys are gene_codes and values
                      are tuples containing BioPython seq_record objects:
                      seq=Seq('????CAGATAAAS'),
                      id='CP13-01_Genus_Species',
                      name='CAD',  # gene_code
                      description='CP13-01',  # voucher code

        ``dataset_str``: output dataset to pass to users.

    """
    def __init__(self, cleaned_data):
        try:
            self.degen_translations = cleaned_data['degen_translations']
        except KeyError:
            self.degen_translations = None

        try:
            self.translations = cleaned_data['translations']
        except KeyError:
            self.translations = None

        self.errors = []
        self.seq_objs = []
        self.minimum_number_of_genes = cleaned_data['number_genes']
        self.aminoacids = cleaned_data['aminoacids']
        self.codon_positions = clean_positions(cleaned_data['positions'])
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
        if len(self.codon_positions) == 2 and \
                '1st' in self.codon_positions and '3rd' in self.codon_positions:
            self.errors = ['Cannot create dataset for only codon positions 1 and 3.']
            self.dataset_file = None
            return ''
        if len(self.codon_positions) == 2 and \
                '2nd' in self.codon_positions and '3rd' in self.codon_positions:
            self.errors = ['Cannot create dataset for only codon positions 2 and 3.']
            self.dataset_file = None
            return ''

        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()
        if self.file_format == 'MEGA':
            fasta = CreateMEGA(self.codon_positions, self.partition_by_positions,
                               self.seq_objs, self.gene_codes, self.voucher_codes,
                               self.file_format, aminoacids=self.aminoacids,
                               degen_translations=self.degen_translations, translations=self.translations)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.errors += fasta.errors
            self.warnings += fasta.warnings
            self.dataset_file = fasta.dataset_file
            self.aa_dataset_file = fasta.aa_dataset_file
            return fasta_dataset

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
                                self.file_format, degen_translations=self.degen_translations,
                                translations=self.translations, aminoacids=self.aminoacids)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            self.dataset_file = fasta.dataset_file
            return fasta_dataset

        if self.file_format == 'PHY':
            phy = CreatePhylip(self.codon_positions, self.partition_by_positions,
                               self.seq_objs, self.gene_codes, self.voucher_codes,
                               self.file_format, self.outgroup, self.voucher_codes_metadata,
                               self.minimum_number_of_genes, self.aminoacids,
                               degen_translations=self.degen_translations, translations=self.translations)
            phylip_dataset = phy.from_seq_objs_to_dataset()
            self.warnings += phy.warnings
            self.dataset_file = phy.dataset_file
            self.charset_block = phy.charset_block
            return phylip_dataset

        if self.file_format == 'TNT':
            tnt = CreateTNT(self.codon_positions, self.partition_by_positions,
                            self.seq_objs, self.gene_codes, self.voucher_codes,
                            self.file_format, self.outgroup, self.voucher_codes_metadata,
                            self.minimum_number_of_genes, self.aminoacids,
                            degen_translations=self.degen_translations, translations=self.translations)
            tnt_dataset = tnt.from_seq_objs_to_dataset()
            self.warnings += tnt.warnings
            self.dataset_file = tnt.dataset_file
            return tnt_dataset

        if self.file_format == 'NEXUS':
            if self.partition_by_positions == 'ONE':
                self.partition_by_positions = 'by gene'
            elif self.partition_by_positions == 'EACH':
                self.partition_by_positions = 'by codon position'

            try:
                dataset = Dataset(self.seq_objs, format='NEXUS', partitioning=self.partition_by_positions,
                                  codon_positions=self.codon_positions[0])
            except ValueError:
                msg = 'You need to specify the reading frame of all genes to do the partitioning by codon positions'
                self.errors.append(msg)
                return ''

            dataset_handler = DatasetHandler(dataset.dataset_str, self.file_format)
            self.dataset_file = dataset_handler.dataset_file
            print(dataset.dataset_str)
            return dataset.dataset_str

    def create_seq_objs(self):
        """Generate a list of SeqRecord-expanded objects.
        """
        sorted_gene_codes = sorted(list(self.gene_codes), key=str.lower)
        our_taxon_names = self.get_taxon_names_for_taxa()
        all_seqs = self.get_all_sequences()

        for gene_code in sorted_gene_codes:
            for code in self.voucher_codes:
                seq_obj = self.build_seq_obj(code, gene_code, our_taxon_names, all_seqs)
                if seq_obj is None:
                    self.warnings += ['Could not find voucher {0}'.format(code)]
                    continue
                self.seq_objs.append(seq_obj)

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
        return seqs_dict

    def build_seq_obj(self, code, gene_code, our_taxon_names, all_seqs):
        """
        Builds a SeqRecordExpanded object. I cannot be built, returns None.
        """
        this_voucher_seqs = self.extract_sequence_from_all_seqs_in_db(all_seqs, code, gene_code)

        if this_voucher_seqs == '?':
            seq = '?' * self.gene_codes_metadata[gene_code]['length']
        else:
            seq = self.create_seq_record(this_voucher_seqs)

        seq_record = SeqRecordExpanded(seq)

        if code in our_taxon_names:
            seq_record.voucher_code = code
            seq_record.taxonomy = our_taxon_names[code]
            seq_record.gene_code = gene_code
            seq_record.reading_frame = self.gene_codes_metadata[gene_code]['reading_frame']
            seq_record.table = self.gene_codes_metadata[gene_code]['genetic_code']
            return seq_record
        else:
            return None

    def extract_sequence_from_all_seqs_in_db(self, all_seqs, code, gene_code):
        try:
            voucher_sequences = all_seqs[code]
        except KeyError:
            self.warnings += ['Could not find sequences for voucher {} and gene_code {}'.format(code, gene_code)]
            return '?'

        try:
            this_voucher_seqs = voucher_sequences[gene_code]
        except KeyError:
            self.warnings += ['Could not find sequences for voucher {} and gene_code {}'.format(code, gene_code)]
            return '?'
        return this_voucher_seqs

    def create_seq_record(self, s):
        """
        Adds ? if the sequence is not long enough
        :param s:
        :return: str.
        """
        gene_code = s['gene_code']
        length = self.gene_codes_metadata[gene_code]['length']
        sequence = s['sequences']
        length_difference = length - len(sequence)

        sequence += '?' * length_difference
        return sequence

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
        queryset = Genes.objects.all().values('gene_code', 'length', 'reading_frame', 'genetic_code')
        gene_codes_metadata = dict()
        for i in queryset:
            gene_code = i['gene_code']
            gene_codes_metadata[gene_code] = {
                'length': i['length'],
                'reading_frame': i['reading_frame'],
                'genetic_code': i['genetic_code'],
            }
        return gene_codes_metadata
