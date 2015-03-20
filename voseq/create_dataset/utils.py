import collections
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import flatten_taxon_names_dict
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers
from .dataset import Dataset


class CreateFasta(Dataset):
    pass


class CreateTNT(Dataset):
    def __init__(self, *args, **kwargs):
        super(CreateTNT, self).__init__(*args, **kwargs)
        self.number_taxa = len(self.voucher_codes)
        self.number_chars = self.get_number_chars_from_gene_codes()

    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
        on needed dataset.
        """
        out = 'nstates dna;\nxread\n'
        out += str(self.number_chars) + ' ' + str(self.number_taxa)

        for i in partitions:
            out += '\n'
            out += '\n'.join(i)

        out += '\n;\nproc/;'
        return out.strip()

    def get_number_chars_from_gene_codes(self):
        chars = 0

        res = Genes.objects.all().values('gene_code', 'length')
        gene_lengths = {i['gene_code']: i['length'] for i in res}

        for gene in self.gene_codes:
            chars += gene_lengths[gene]
        return chars


class CreateNEXUS(Dataset):
    def __init__(self, *args, **kwargs):
        super(CreateNEXUS, self).__init__(*args, **kwargs)
        self.gene_codes_and_lengths = None
        self.number_taxa = len(self.voucher_codes)
        self.number_chars = None
        self.vouchers_to_drop = None

    def get_charset_block(self):
        charset_block = []

        bp_count_start = 0
        bp_count_end = 0
        self.gene_codes.sort()
        for gene in self.gene_codes_and_lengths:
            bp_count_end += self.gene_codes_and_lengths[gene]
            line = '    charset ' + gene + ' = ' + str(
                bp_count_start + 1) + '-' + str(bp_count_end) + ';'
            bp_count_start += self.gene_codes_and_lengths[gene]
            charset_block.append(line)
        return charset_block

    def get_partitions_block(self):
        line = 'partition GENES = ' + str(len(self.gene_codes_and_lengths))
        line += ': ' + ', '.join([i for i in self.gene_codes_and_lengths]) + ';\n'
        line += '\nset partition = GENES;\n'
        return [line]

    def get_final_block(self):
        block = "set autoclose=yes;"
        if self.outgroup != '':
            block += '\noutgroup ' + self.voucher_codes_metadata[self.outgroup] + ";"
        block += """
prset applyto=(all) ratepr=variable brlensp=unconstrained:Exp(100.0) shapepr=exp(1.0) tratiopr=beta(2.0,1.0);
lset applyto=(all) nst=mixed rates=gamma [invgamma];
unlink statefreq=(all);
unlink shape=(all) revmat=(all) tratio=(all) [pinvar=(all)];
mcmc ngen=10000000 printfreq=1000 samplefreq=1000 nchains=4 nruns=2 savebrlens=yes [temp=0.11];
 sump relburnin=yes [no] burninfrac=0.25 [2500];
 sumt relburnin=yes [no] burninfrac=0.25 [2500] contype=halfcompat [allcompat];
END;
"""
        return [block.strip()]

    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        out = [
            '#NEXUS\n',
            'BEGIN DATA;',
            'DIMENSIONS NTAX=' + str(self.number_taxa - len(self.vouchers_to_drop)) + ' NCHAR=' + str(self.number_chars) + ';',
            'FORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;',
            'MATRIX',
        ]

        for partition in partitions:
            for i in partition:
                if i.split(' ')[0] not in self.vouchers_to_drop:
                    out += [i]

        out += [';\nEND;']
        out += ['\nbegin mrbayes;']
        out += self.get_charset_block()
        out += self.get_partitions_block()
        out += self.get_final_block()
        return '\n'.join(out)

    def get_number_chars_from_partition_list(self, partitions):
        chars = 0

        gene_codes_and_lengths = collections.OrderedDict()

        gene_code = ''
        for item in partitions[0]:
            if item.startswith('\n'):
                gene_code = item.strip().replace('[', '').replace(']', '')
                continue
            if gene_code != '':
                first_entry = re.sub('\s+', ' ', item)
                voucher, sequence = first_entry.split(' ')
                chars += len(sequence)
                gene_codes_and_lengths[gene_code] = len(sequence)
                gene_code = ''
        self.gene_codes_and_lengths = gene_codes_and_lengths
        self.number_chars = chars


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
        self.codon_positions = cleaned_data['positions']
        self.file_format = cleaned_data['file_format']
        self.partition_by_positions = cleaned_data['partition_by_positions']
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.taxon_names = cleaned_data['taxon_names']
        self.voucher_codes_metadata = dict()
        self.gene_codes_metadata = self.get_gene_codes_metadata()
        self.warnings = []
        self.outgroup = cleaned_data['outgroup']
        self.dataset_str = self.create_dataset()

    def create_dataset(self):
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

        if self.file_format == 'FASTA':
            fasta = CreateFasta(self.codon_positions, self.partition_by_positions,
                                self.seq_objs, self.gene_codes, self.voucher_codes,
                                self.file_format)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            return fasta_dataset

        if self.file_format == 'TNT':
            tnt = CreateTNT(self.codon_positions, self.partition_by_positions,
                            self.seq_objs, self.gene_codes, self.voucher_codes, self.file_format, self.outgroup)
            tnt_dataset = tnt.from_seq_objs_to_dataset()
            self.warnings += tnt.warnings
            return tnt_dataset

        if self.file_format == 'NEXUS':
            nexus = CreateNEXUS(self.codon_positions, self.partition_by_positions,
                                self.seq_objs, self.gene_codes, self.voucher_codes,
                                self.file_format, self.outgroup, self.voucher_codes_metadata,
                                self.minimum_number_of_genes)
            nexus_dataset = nexus.from_seq_objs_to_dataset()
            self.warnings += nexus.warnings
            return nexus_dataset

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

    def create_seq_objs(self):
        """Generate a list of sequence objects. Also takes into account the
        genes passed as geneset.

        Returns:
            list of sequence objects as produced by BioPython.

        """
        # We might need to update our list of vouches and genes
        vouchers_found = set()
        gene_codes = set()

        our_taxon_names = self.get_taxon_names_for_taxa()

        all_seqs = Sequences.objects.all().values('code_id', 'gene_code', 'sequences').order_by('code_id')
        for s in all_seqs:
            code = s['code_id']
            gene_code = s['gene_code']
            if code in self.voucher_codes and gene_code in self.gene_codes:
                vouchers_found.add(code)
                gene_codes.add(gene_code)

                seq_obj = self.create_seq_record(s)
                seq_obj.id = flatten_taxon_names_dict(our_taxon_names[code])
                if 'GENECODE' in self.taxon_names:
                    seq_obj.id += '_' + gene_code
                seq_obj.name = gene_code
                seq_obj.description = code
                self.voucher_codes_metadata[code] = seq_obj.id

                if gene_code not in self.seq_objs:
                    self.seq_objs[gene_code] = []
                self.seq_objs[gene_code].append(seq_obj)

        vouchers_not_found = set(self.voucher_codes) - vouchers_found
        self.warnings += ['Could not found sequences for voucher %s' % i for i in vouchers_not_found]
        self.voucher_codes = list(vouchers_found)
        self.gene_codes = list(gene_codes)
        self.add_missing_seqs()

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

    def add_missing_seqs(self):
        """
        Loops over the created seq_objects and adds sequences full of ? if
        those where not found in our database.

        Uses the updated lists of voucher_codes and gene_codes
        """
        for gene_code in self.seq_objs:
            for code in self.voucher_codes:
                found = False
                for seq_obj in self.seq_objs[gene_code]:
                    if code == seq_obj.description:
                        found = True

                if found is False:
                    seq = Seq('?' * self.gene_codes_metadata[gene_code])
                    empty_seq_obj = SeqRecord(seq)
                    empty_seq_obj.id = self.voucher_codes_metadata[code]
                    empty_seq_obj.name = gene_code
                    empty_seq_obj.description = code
                    self.seq_objs[gene_code].append(empty_seq_obj)

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
