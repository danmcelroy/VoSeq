from seqrecord_expanded import SeqRecordExpanded
from seqrecord_expanded.exceptions import MissingParameterError
from seqrecord_expanded.exceptions import TranslationErrorMixedGappedSeq
from dataset_creator import Dataset
from Bio.Nexus.Nexus import NexusError

from core import exceptions
from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import clean_positions
from .nexus import DatasetHandler
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers


class CreateDataset(object):
    """Accepts form input to create a dataset in several formats.

    Also codon positions, for list of codes and genes. Takes into account the
    vouchers passed as taxonset.

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
        self.cleaned_data = cleaned_data
        self.translations = None
        self.degen_translations = None
        self.clean_translations()

        self.errors = []
        self.seq_objs = []
        self.minimum_number_of_genes = cleaned_data['number_genes']
        self.aminoacids = cleaned_data['aminoacids']

        try:
            self.codon_positions = clean_positions(cleaned_data['positions'])
        except exceptions.InadequateCodonPositions as e:
            self.codon_positions = None
            self.errors = [e]

        self.file_format = cleaned_data['file_format']
        self.partition_by_positions = cleaned_data['partition_by_positions']
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

    def clean_translations(self):
        if self.cleaned_data['translations']:
            self.degen_translations = self.cleaned_data['degen_translations']
        else:
            # No need to do degen translation
            self.degen_translations = None

    def create_dataset(self):
        if not self.codon_positions:
            return ''

        if self.degen_translations is not None and self.codon_positions != ['ALL']:
            msg = 'Cannot degenerate codons if you have not selected all codon positions'
            self.errors.append(msg)
            return ''
        elif self.degen_translations is not None and self.partition_by_positions != 'by gene':
            msg = 'Cannot degenerate codons if they go to different partitions'
            self.errors.append(msg)
            return ''

        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

        if self.file_format in ['NEXUS', 'GenBankFASTA', 'FASTA', 'MEGA', 'TNT', 'PHYLIP']:
            try:
                dataset = Dataset(self.seq_objs, format=self.file_format, partitioning=self.partition_by_positions,
                                  codon_positions=self.codon_positions[0], aminoacids=self.aminoacids,
                                  degenerate=self.degen_translations, outgroup=self.outgroup)
            except MissingParameterError as e:
                self.errors.append(e)
                return ''
            except ValueError as e:
                self.errors.append(e)
                return ''
            except TranslationErrorMixedGappedSeq as e:
                self.errors.append(e)
                return ''
            except NexusError as e:
                self.errors.append(e.__str__())
                return ''

            self.warnings += dataset.warnings
            dataset_handler = DatasetHandler(dataset.dataset_str, self.file_format)
            self.dataset_file = dataset_handler.dataset_file

            if self.file_format == 'PHYLIP':
                self.charset_block = dataset.extra_dataset_str

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
        """Return sequences as dict of lists containing sequence and related data.

        """
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
        """Builds a SeqRecordExpanded object. I cannot be built, returns None.

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
            self.warnings += ['Could not find sequences for voucher {0} and gene_code {1}'.format(code, gene_code)]
            return '?'

        try:
            this_voucher_seqs = voucher_sequences[gene_code]
        except KeyError:
            self.warnings += ['Could not find sequences for voucher {0} and gene_code {1}'.format(code, gene_code)]
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
