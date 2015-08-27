import collections
from collections import namedtuple
import os
import uuid
import re

from core import utils
from public_interface.models import Genes


class Dataset(object):
    """
    Base class to create datasets from Seq objects into FASTA, TNT formats.
    """
    def __init__(self, codon_positions, partition_by_positions, seq_objs, gene_codes,
                 voucher_codes, file_format, outgroup=None, voucher_codes_metadata=None,
                 minimum_number_of_genes=None, aminoacids=None, degen_translations=None,
                 translations=None):
        self.degen_translations = degen_translations
        self.translations = translations
        self.minimum_number_of_genes = minimum_number_of_genes
        self.outgroup = outgroup
        self.file_format = file_format
        self.codon_positions = codon_positions
        self.aminoacids = aminoacids
        self.partition_by_positions = partition_by_positions
        self.seq_objs = seq_objs
        self.gene_codes = gene_codes
        self.gene_code_descriptions = self.get_gene_code_descriptions()
        self.voucher_codes = voucher_codes
        self.vouchers_to_drop = None
        self.number_taxa = len(self.voucher_codes)
        self.reading_frames = self.get_reading_frames()
        self.voucher_codes_metadata = voucher_codes_metadata
        self.warnings = []
        self.partition_list = None

        self.cwd = os.path.dirname(__file__)
        self.guid = self.make_guid()
        self.dataset_file = os.path.join(self.cwd,
                                         'dataset_files',
                                         self.file_format + '_' + self.guid + '.txt',
                                         )
        self.aa_dataset_file = os.path.join(self.cwd,
                                            'dataset_files',
                                            self.file_format + '_aa_' + self.guid + '.txt',
                                            )

    def save_dataset_to_file(self, dataset_str):
        with open(self.dataset_file, 'w') as handle:
            handle.write(dataset_str)

    def save_aa_dataset_to_file(self, aa_dataset_str):
        with open(self.aa_dataset_file, 'w') as handle:
            handle.write(aa_dataset_str)

    def make_guid(self):
        return uuid.uuid4().hex

    def get_gene_for_current_partition(self, gene_models, out,
                                       partitions_incorporated, voucher_code):
        ThisGeneAndPartition = namedtuple('ThisGeneAndPartition', ['this_gene',
                                                                   'this_gene_model',
                                                                   'partitions_incorporated',
                                                                   'out'])
        gene = voucher_code.replace('[', '').replace(']', '').strip()
        ThisGeneAndPartition.this_gene = gene
        ThisGeneAndPartition.this_gene_model = self.get_gene_model_from_gene_id(ThisGeneAndPartition.this_gene, gene_models)
        ThisGeneAndPartition.partitions_incorporated = partitions_incorporated + 1
        if self.file_format == 'NEXUS':
            ThisGeneAndPartition.out = out + ['\n[' + gene + ']']
        else:
            ThisGeneAndPartition.out = out + ['\n']
        return ThisGeneAndPartition

    def get_gene_model_from_gene_id(self, this_gene, gene_models):
        for i in gene_models:
            if i['gene_code'] == this_gene:
                return i

    def get_number_chars_from_partition_list(self, partitions):
        chars = 0
        gene_codes_and_lengths = collections.OrderedDict()

        i = 0
        gene_code = ''
        for item in partitions[0]:
            if item.startswith('\n'):
                if self.file_format == 'NEXUS' or self.file_format == 'PHY':
                    gene_code = item.strip().replace('[', '').replace(']', '')
                    continue
                if self.file_format == 'TNT':
                    gene_code = 'dummy' + str(i)
                    i += 1
                    continue
            if gene_code != '':
                first_entry = re.sub('\s+', ' ', item)
                sequence = first_entry.split(' ')[1]
                chars += len(sequence)
                gene_codes_and_lengths[gene_code] = len(sequence)
                gene_code = ''
        self.gene_codes_and_lengths = gene_codes_and_lengths
        self.number_chars = chars

    def get_number_of_genes_for_taxa(self, partitions):
        number_of_genes_for_taxa = dict()
        vouchers_to_drop = set()

        gene_code = ''
        for item in partitions[0]:
            if item.startswith('\n'):
                if self.file_format == 'NEXUS' or self.file_format == 'PHY':
                    gene_code = item.strip().replace('[', '').replace(']', '')
                    continue
                if self.file_format == 'TNT':
                    gene_code = 'dummy'
                    continue
            if gene_code != '':
                entry = re.sub('\s+', ' ', item)
                voucher, sequence = entry.split(' ')

                if voucher not in number_of_genes_for_taxa:
                    number_of_genes_for_taxa[voucher] = 0

                sequence = sequence.replace('?', '')
                if sequence != '':
                    number_of_genes_for_taxa[voucher] += 1

        if self.minimum_number_of_genes is None:
            self.vouchers_to_drop = []
        else:
            for voucher in number_of_genes_for_taxa:
                if number_of_genes_for_taxa[voucher] < self.minimum_number_of_genes:
                    vouchers_to_drop.add(voucher)
            self.vouchers_to_drop = vouchers_to_drop

    def get_gene_code_descriptions(self):
        gene_metadata = dict()
        genes = Genes.objects.all().values('gene_code', 'description')
        for i in genes:
            gene_code = i['gene_code']
            gene_metadata[gene_code] = i['description']
        return gene_metadata

    def get_reading_frames(self):
        """
        :return: dict of gene_code: reading_frame
        """
        reading_frames = dict()
        genes = Genes.objects.all().values('gene_code', 'reading_frame')
        for gene in genes:
            gene_code = gene['gene_code']
            if gene_code in self.gene_codes:
                reading_frames[gene_code] = gene['reading_frame']
        return reading_frames

    def split_sequence_in_codon_positions(self, gene_code, seq):
        """
        Retuns tuple of nucleotides based on codon positions.

        :param gene_code: as lower case
        :param seq: as BioPython seq object.
        :return: tuple of seq strings.
        """
        reading_frame = int(self.reading_frames[gene_code])

        # This is the BioPython way to get codon positions
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc19
        if reading_frame == 1:
            first_position = seq[::3]
            second_position = seq[1::3]
            third_position = seq[2::3]
        elif reading_frame == 2:
            first_position = seq[1::3]
            second_position = seq[2::3]
            third_position = seq[::3]
        elif reading_frame == 3:
            first_position = seq[2::3]
            second_position = seq[::3]
            third_position = seq[1::3]
        return first_position, second_position, third_position

    def convert_lists_to_dataset(self, partitions):
        """
        Method to override in order to add headers and footers depending
        on needed dataset.
        """
        out = ''
        for partition in partitions:
            this_gene = ''
            for line in partition:
                line_contents = line.split('\n')
                if line_contents[1] == '--------------------':
                    this_gene_partition = line_contents[0].replace('>', '').strip()
                    this_gene = Genes.objects.filter(gene_code=this_gene_partition).values()
                    try:
                        this_gene = this_gene[0]
                    except IndexError:
                        pass
                    out += '\n'
                    out += line
                else:
                    line_contents = line.split('\n')
                    taxon = line_contents[0]
                    voucher_code = taxon.replace('>', '').split('_')[0]
                    sequence = line_contents[1]

                    if self.aminoacids is not True and this_gene != '' and \
                            self.translations is True:
                        sequence = self.degenerate(sequence, this_gene)

                    if self.aminoacids is True:
                        sequence = self.translate_this_sequence(
                            sequence,
                            this_gene,
                            voucher_code,
                        )

                    out += '\n'
                    out += '\n'.join([taxon, sequence])

        dataset_str = out.strip()
        self.save_dataset_to_file(dataset_str)
        return dataset_str

    def degenerate(self, seq, gene_model):
        if self.translations is None or self.translations is False:
            return seq

        if self.degen_translations in ['NORMAL', 'S', 'Z', 'SZ']:
            if self.partition_by_positions != 'ONE':
                self.warnings.append(
                    'Cannot degenerate codons if they go to different partitions.'
                )
                return ''

            if 'ALL' in self.codon_positions or (
                    '1st' in self.codon_positions and
                    '2nd' in self.codon_positions and
                    '3rd' in self.codon_positions):
                degenerated = utils._degenerate(gene_model, seq, self.degen_translations)
                return degenerated
            else:
                self.warnings.append(
                    'Cannot degenerate codons if they you have not selected all codon positions.'
                )
                return ''

    def translate_this_sequence(self, sequence, this_gene_model, voucher_code):
        aa_sequence = ''
        if this_gene_model['genetic_code'] is None or this_gene_model['reading_frame'] is None:
            self.warnings.append(
                "Cannot translate gene %s sequences into aminoacids."
                " You need to define reading_frame and/or genetic_code." %
                this_gene_model['gene_code'])
        else:
            aa_sequence, warning = utils.translate_to_protein(this_gene_model,
                                                              sequence, '',
                                                              voucher_code,
                                                              self.file_format)
            if warning != '':
                self.warnings.append(warning)
        sequence = aa_sequence

        if sequence == '':
            sequence = '?'
        return sequence

    def get_codons_in_each_partition(self, codons):
        partition_list = ()
        codon_descriptions = []
        codon_pos = []
        for i in codons:
            partition_list += ([],)
            if i == '1st':
                codon_descriptions.append('_1st_codon')
                codon_pos.append(0)
            if i == '2nd':
                codon_descriptions.append('_2nd_codon')
                codon_pos.append(1)
            if i == '3rd':
                codon_descriptions.append('_3rd_codon')
                codon_pos.append(2)
        for gene_code in self.seq_objs:
            this_gene = None

            for seq_record in self.seq_objs[gene_code]:
                if self.reading_frames[gene_code] is None:
                    self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                         "it cannot be included in your dataset." % gene_code)
                    continue

                if this_gene is None:
                    this_gene = seq_record.name

                    for i in range(len(codon_descriptions)):
                        seq_str = self.get_gene_divisor(this_gene, codon_descriptions[i])
                        partition_list[i].append(seq_str)

                codons = self.split_sequence_in_codon_positions(this_gene,
                                                                seq_record.seq)
                for i in range(len(codon_pos)):
                    seq_str = self.format_record_id_and_seq_for_dataset(seq_record.id, codons[codon_pos[i]])
                    partition_list[i].append(seq_str)
        return partition_list

    def get_gene_divisor(self, this_gene, codon_description=None):
        if self.file_format == 'FASTA':
            seq_str = '>' + this_gene
            if codon_description is not None:
                seq_str += codon_description
            seq_str += '\n' + '--------------------'

        if self.file_format == 'MEGA':
            seq_str = '\n[%s]' % this_gene

        if self.file_format == 'GenbankFASTA':
            seq_str = '\n[%s]' % this_gene

        if self.file_format == 'PHY':
            seq_str = '\n[%s]' % this_gene

        if self.file_format == 'TNT':
            seq_str = '\n[%s]' % this_gene

        if self.file_format == 'NEXUS':
            seq_str = '\n[%s]' % this_gene
        return seq_str

    def format_record_id_and_seq_for_dataset(self, seq_record_id, seq_record_seq):
        if self.file_format == 'GenbankFASTA':
            seq_record_id = seq_record_id.split(' ')
            code = seq_record_id[0]
            gene = seq_record_id[1]
            taxon = seq_record_id[2].replace('_', ' ')
            gene_description = self.gene_code_descriptions[gene]

            seq_str = '>' + code + '_' + gene + ' ' + '[organism=' + taxon + '] [specimen-voucher=' + code + '] ' + taxon
            seq_str += ' ' + gene_description + '\n' + str(seq_record_seq)
        if self.file_format == 'FASTA':
            seq_str = '>' + seq_record_id + '\n' + str(seq_record_seq)
        if self.file_format == 'NEXUS' or \
                self.file_format == 'PHY' or \
                self.file_format == 'MEGA' or \
                self.file_format == 'TNT':
            seq_str = str(seq_record_id).ljust(55) + str(seq_record_seq)
        return seq_str

    def get_codons_in_one_partition(self, codons):
        partition_list = ([],)
        codon_pos = []
        for i in codons:
            if i == '1st':
                codon_pos.append(0)
            if i == '2nd':
                codon_pos.append(1)
            if i == '3rd':
                codon_pos.append(2)
        for gene_code in self.seq_objs:
            this_gene = None

            for seq_record in self.seq_objs[gene_code]:
                if self.reading_frames[gene_code] is None:
                    self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                         "it cannot be included in your dataset." % gene_code)
                    continue

                if this_gene is None:
                    this_gene = seq_record.name

                    gene_divisor = self.get_gene_divisor(this_gene)
                    partition_list[0].append(gene_divisor)

                codons = self.split_sequence_in_codon_positions(this_gene,
                                                                seq_record.seq)

                codon_seqs = str(utils.chain_and_flatten([codons[i] for i in codon_pos]))
                seq_str = self.format_record_id_and_seq_for_dataset(seq_record.id, codon_seqs)

                partition_list[0].append(seq_str)
        self.partition_list = partition_list
        return partition_list

    def from_seq_objs_to_dataset(self):
        """Take a list of BioPython's sequence objects and return a FASTA string

        Returns:
            The FASTA string will contain the gene_code at the top as if it were
            another FASTA gene sequence.

        """
        if '1st' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'EACH' in self.partition_by_positions and \
                self.file_format == 'PHY':
            self.partition_list = self.get_codons_in_one_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'EACH' in self.partition_by_positions and \
                self.file_format != 'PHY':
            self.partition_list = self.get_codons_in_each_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_one_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if ('ALL' in self.codon_positions or
                ('1st' in self.codon_positions and '2nd' in self.codon_positions and '3rd' in self.codon_positions)) \
                and 'EACH' in self.partition_by_positions and self.file_format in ['PHY', 'NEXUS']:
            self.get_all_codon_positions_in_one_partition()
            return self.convert_lists_to_dataset(self.partition_list)

        if ('ALL' in self.codon_positions or
                ('1st' in self.codon_positions and '2nd' in self.codon_positions and '3rd' in self.codon_positions)) \
                and 'EACH' in self.partition_by_positions and self.file_format != 'PHY':
            self.partition_list = self.get_codons_in_each_partition(['1st', '2nd', '3rd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if ('ALL' in self.codon_positions or
                ('1st' in self.codon_positions and '2nd' in self.codon_positions and '3rd' in self.codon_positions)) \
                and 'ONE' in self.partition_by_positions:
            self.get_all_codon_positions_in_one_partition()
            return self.convert_lists_to_dataset(self.partition_list)

        if 'ALL' in self.codon_positions and '1st2nd_3rd' in self.partition_by_positions and \
                self.file_format in ['PHY', 'NEXUS']:
            self.get_all_codon_positions_in_one_partition()
            return self.convert_lists_to_dataset(self.partition_list)

        if 'ALL' in self.codon_positions and '1st2nd_3rd' in self.partition_by_positions and self.file_format != 'PHY':
            self.partition_list = ([], [],)

            for gene_code in self.seq_objs:
                if self.reading_frames[gene_code] is None:
                    self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                         "it cannot be included in your dataset." % gene_code)
                    continue

                this_gene = None
                for seq_record in self.seq_objs[gene_code]:
                    if this_gene is None:
                        this_gene = seq_record.name

                        seq_str = '>' + this_gene + '_1st_2nd_codons\n' + '--------------------'
                        self.partition_list[0].append(seq_str)
                        seq_str = '>' + this_gene + '_3rd_codon\n' + '--------------------'
                        self.partition_list[1].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(utils.chain_and_flatten([codons[0], codons[1]]))
                    self.partition_list[0].append(seq_str)

                    seq_str = '>' + seq_record.id + '\n' + str(codons[2])
                    self.partition_list[1].append(seq_str)
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            self.partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions and \
                self.file_format == 'PHY':
            self.partition_list = self.get_codons_in_one_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions and \
                self.file_format != 'PHY':
            self.partition_list = ([],)
            for gene_code in self.seq_objs:
                this_gene = None
                for seq_record in self.seq_objs[gene_code]:
                    if self.reading_frames[gene_code] is None:
                        self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                             "it cannot be included in your dataset." % gene_code)
                        continue

                    if this_gene is None:
                        this_gene = seq_record.name

                        seq_str = '>' + this_gene + '_1st_2nd_codons\n' + '--------------------'
                        self.partition_list[0].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(utils.chain_and_flatten([codons[0], codons[1]]))
                    self.partition_list[0].append(seq_str)
            return self.convert_lists_to_dataset(self.partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            self.partition_list = ([], [],)
            for gene_code in self.seq_objs:
                this_gene = None
                for seq_record in self.seq_objs[gene_code]:
                    if self.reading_frames[gene_code] is None:
                        self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                             "it cannot be included in your dataset." % gene_code)
                        continue

                    if this_gene is None:
                        this_gene = seq_record.name

                        seq_str = '>' + this_gene + '_1st_2nd_codons\n' + '--------------------'
                        self.partition_list[0].append(seq_str)

                        seq_str = '>' + this_gene + '_3rd_codon\n' + '--------------------'
                        self.partition_list[1].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(utils.chain_and_flatten([codons[0], codons[1]]))
                    self.partition_list[0].append(seq_str)

                    seq_str = '>' + seq_record.id + '\n' + str(codons[2])
                    self.partition_list[1].append(seq_str)
            return self.convert_lists_to_dataset(self.partition_list)

    def get_all_codon_positions_in_one_partition(self):
        self.partition_list = ([],)
        for gene_code in self.seq_objs:
            this_gene = None
            for seq_record in self.seq_objs[gene_code]:
                if this_gene is None:
                    this_gene = seq_record.name

                    seq_str = self.get_gene_divisor(this_gene)
                    self.partition_list[0].append(seq_str)

                if self.file_format == 'GenbankFASTA':
                    seq_str = self.format_record_id_and_seq_for_dataset(
                        seq_record.description + ' ' +
                        this_gene + ' ' +
                        seq_record.id.replace(seq_record.description + '_',
                                              ''), seq_record.seq)
                else:
                    seq_str = self.format_record_id_and_seq_for_dataset(
                        seq_record.id, seq_record.seq)
                self.partition_list[0].append(seq_str)

    def make_charset_block(self, gene_codes_and_lengths):
        gene_codes_list = sorted(list(gene_codes_and_lengths), key=str.lower)
        charset_block = self.generate_charset_block(gene_codes_and_lengths, gene_codes_list)
        self.charset_block = "\n".join(charset_block)

    def generate_charset_block(self, gene_codes_and_lengths, gene_codes_list):
        """Basic charset block. All codon positions and one partition for gene.
        """
        bp_count_start = 0
        bp_count_end = 0
        charset_block = []
        for gene in gene_codes_list:
            bp_count_end += gene_codes_and_lengths[gene]

            line = self.make_charset_line(bp_count_start, bp_count_end, gene)

            bp_count_start += gene_codes_and_lengths[gene]
            charset_block.append(line)
        return charset_block

    def make_charset_line(self, count_start, count_end, gene):
        """To be used when making the charset block. It writes one line per partition.
        """
        if self.file_format == 'PHY':
            prefix = 'DNA,'
            suffix = ''
        elif self.file_format == 'NEXUS':
            prefix = '    charset'
            suffix = ';'

        if self.partition_by_positions == 'ONE':
            if 'ALL' in self.codon_positions:
                line = '{} {} = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = '{} {}_pos1 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = '{} {}_pos2 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = '{} {}_pos3 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 2 and \
                    '1st' in self.codon_positions and \
                    '2nd' in self.codon_positions:
                line = '{} {}_pos12 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line

        elif self.partition_by_positions == '1st2nd_3rd':
            if 'ALL' in self.codon_positions:
                line = ''
                if self.reading_frames[gene] == 1:
                    line = '{} {}_pos12 = '.format(prefix, gene)
                    line += '{}-{}\\3'.format(count_start + 1, count_end)
                    line += ', {}-{}\\3{}'.format(count_start + 2, count_end, suffix)
                    line += '\n{} {}_pos3 = '.format(prefix, gene)
                    line += '{}-{}\\3{}'.format(count_start + 3, count_end, suffix)
                elif self.reading_frames[gene] == 2:
                    line = '{} {}_pos12 = '.format(prefix, gene)
                    line += '{}-{}\\3'.format(count_start + 2, count_end)
                    line += ', {}-{}\\3{}'.format(count_start + 3, count_end, suffix)
                    line += '\n{} {}_pos3 = '.format(prefix, gene)
                    line += '{}-{}\\3{}'.format(count_start + 1, count_end, suffix)
                elif self.reading_frames[gene] == 3:
                    line = '{} {}_pos12 = '.format(prefix, gene)
                    line += '{}-{}\\3'.format(count_start + 3, count_end)
                    line += ', {}-{}\\3{}'.format(count_start + 1, count_end, suffix)
                    line += '\n{} {}_pos3 = '.format(prefix, gene)
                    line += '{}-{}\\3{}'.format(count_start + 2, count_end, suffix)
                return line

            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = '{} {}_pos1 = {}-{}'.format(prefix, gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = '{} {}_pos2 = {}-{}'.format(prefix, gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = '{} {}_pos3 = {}-{}'.format(prefix, gene, count_start + 1, count_end)
                return line

            if len(self.codon_positions) == 2 and \
                    '1st' in self.codon_positions and \
                    '2nd' in self.codon_positions:
                line = '{} {}_pos12 = {}-{}'.format(prefix, gene, count_start + 1, count_end)
                return line

        elif self.partition_by_positions == 'EACH':
            if 'ALL' in self.codon_positions:
                line = ''
                if self.reading_frames[gene] == 1:
                    line = '{} {}_pos1 = {}-{}\\3{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                    line += '\n{} {}_pos2 = {}-{}\\3{}'.format(prefix, gene, count_start + 2, count_end, suffix)
                    line += '\n{} {}_pos3 = {}-{}\\3{}'.format(prefix, gene, count_start + 3, count_end, suffix)
                elif self.reading_frames[gene] == 2:
                    line = '{} {}_pos1 = {}-{}\\3{}'.format(prefix, gene, count_start + 2, count_end, suffix)
                    line += '\n{} {}_pos2 = {}-{}\\3{}'.format(prefix, gene, count_start + 3, count_end, suffix)
                    line += '\n{} {}_pos3 = {}-{}\\3{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                elif self.reading_frames[gene] == 3:
                    line = '{} {}_pos1 = {}-{}\\3{}'.format(prefix, gene, count_start + 3, count_end, suffix)
                    line += '\n{} {}_pos2 = {}-{}\\3{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                    line += '\n{} {}_pos3 = {}-{}\\3{}'.format(prefix, gene, count_start + 2, count_end, suffix)
                return line

            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = '{} {}_pos1 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = '{} {}_pos2 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = '{} {}_pos3 = {}-{}{}'.format(prefix, gene, count_start + 1, count_end, suffix)
                return line

            if len(self.codon_positions) == 2 and \
                    '1st' in self.codon_positions and \
                    '2nd' in self.codon_positions:
                line = 'DNA, {}_pos1 = {}-{}\\2{}'.format(gene, count_start + 1, count_end, suffix)
                line += '\nDNA, {}_pos2 = {}-{}\\2{}'.format(gene, count_start + 2, count_end, suffix)
                return line


class CreateFasta(Dataset):
    pass


class CreateMEGA(Dataset):
    def convert_lists_to_dataset(self, partitions):
        sequence_dict = dict()
        for partition in partitions:
            this_gene = ''
            for line in partition:
                if line.startswith('\n['):
                    this_gene_partition = line.replace('[', '').replace(']', '').strip()
                    this_gene = Genes.objects.filter(gene_code=this_gene_partition).values()
                    try:
                        this_gene = this_gene[0]
                    except KeyError:
                        pass
                elif line.startswith('>'):
                    voucher_code = line.splitlines()[0].replace('>', '')
                    sequence = line.splitlines()[1]
                    line = '{} {}'.format(voucher_code, sequence)

                line = line.split(' ')
                taxon = line[0].replace('?', '')
                voucher_code = taxon.split('_')[0]
                sequence = line[-1]

                if self.partition_by_positions != 'ONE' and self.translations is True:
                    self.warnings.append('Cannot degenerate codons if they go to different partitions.')
                    continue

                if self.aminoacids is True:
                    sequence = self.translate_this_sequence(
                        sequence,
                        this_gene,
                        voucher_code,
                    )

                if self.aminoacids is not True and this_gene != '' and \
                        self.translations is True:
                    sequence = self.degenerate(sequence, this_gene)

                if taxon not in sequence_dict:
                    sequence_dict[taxon] = ''
                    sequence_dict[taxon] += sequence
                else:
                    sequence_dict[taxon] += sequence

        out = ''
        for k, v in sequence_dict.items():
            out += '\n#' + k + '\n' + v
        dataset_str = '#MEGA\n!TITLE title;\n\n' + out.strip()
        self.save_dataset_to_file(dataset_str)
        return dataset_str


class CreateGenbankFasta(Dataset):
    def convert_lists_to_dataset(self, partitions):
        """
        Overridden method from base class in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        out = []
        aa_out = []

        gene_models = Genes.objects.all().values()

        partitions_incorporated = 0
        for partition in partitions:
            for i in partition:
                i = i.strip()
                if i.startswith('['):
                    this_gene = i.replace('[', '').replace(']', '').strip()
                    this_gene_model = self.get_gene_model_from_gene_id(this_gene, gene_models)
                    partitions_incorporated += 1
                    out += ['\n']
                elif i.startswith('>'):
                    try:
                        voucher_code = re.search('specimen-voucher=(.+)]', i).groups()[0]
                    except AttributeError:
                        continue

                    if voucher_code not in self.vouchers_to_drop:
                        line = i.split('\n')
                        if len(line) > 1:
                            sequence = line[-1]

                            if this_gene_model['genetic_code'] is None or this_gene_model['reading_frame'] is None:
                                self.warnings.append("Cannot translate gene %s sequences into aminoacids."
                                                     " You need to define reading_frame and/or genetic_code." % this_gene_model['gene_code'])
                                continue
                            else:
                                aa_sequence, warning = utils.translate_to_protein(this_gene_model, sequence, '', voucher_code, self.file_format)
                                if warning != '':
                                    self.warnings.append(warning)

                            if aa_sequence.strip() == '':
                                self.warnings.append("Sequence for %s %s was empty" % (voucher_code, this_gene))

                            sequence = utils.strip_question_marks(sequence)[0]
                            out += [line[0] + '\n' + sequence.replace('?', 'N') + '\n']
                            aa_out += [line[0] + '\n' + aa_sequence + '\n']

        dataset_str = ''.join(out)
        aa_dataset_str = ''.join(aa_out)
        self.save_dataset_to_file(dataset_str)
        self.save_aa_dataset_to_file(aa_dataset_str)
        return dataset_str, aa_dataset_str


class CreateTNT(Dataset):
    def convert_lists_to_dataset(self, partitions):
        """
        Overridden method from base class in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        gene_models = Genes.objects.all().values()
        gene_codes_and_lengths = collections.OrderedDict()

        out = []

        outgroup_sequences = []
        for partition in partitions:
            for i in partition:
                voucher = i.split(' ')[0]
                if voucher not in self.vouchers_to_drop:
                    if self.outgroup in voucher:
                        outgroup_sequences.append(i)
                        continue

        partitions_incorporated = 0
        partition_count = 0
        for partition in partitions:
            for i in partition:
                voucher_code = i.split(' ')[0]
                if voucher_code.startswith('\n'):
                    ThisGeneAndPartition = self.get_gene_for_current_partition(
                        gene_models, out, partitions_incorporated, voucher_code
                    )
                    out += ['\n[&dna]']
                    if self.outgroup != '':
                        line = outgroup_sequences[partition_count]
                        tmp_line = line.split(' ')
                        sequence = tmp_line[-1]

                        if self.aminoacids is True:
                            sequence = self.translate_this_sequence(
                                sequence,
                                ThisGeneAndPartition.this_gene_model,
                                voucher_code
                            )
                            out += [tmp_line[0].ljust(55, ' ') + sequence]
                        else:
                            out += [outgroup_sequences[partition_count]]
                        partition_count += 1
                elif voucher_code.startswith('>'):
                    voucher_code = i.splitlines()[0].replace('>', '')
                    sequence = i.splitlines()[1]
                    i = '{} {}'.format(voucher_code, sequence)
                    ThisGeneAndPartition = self.get_gene_for_current_partition(
                        gene_models, out, partitions_incorporated, voucher_code
                    )

                if voucher_code not in self.vouchers_to_drop:
                    line = i.split(' ')
                    if len(line) > 1:
                        sequence = line[-1]

                        if self.aminoacids is True:
                            sequence = self.translate_this_sequence(
                                sequence,
                                ThisGeneAndPartition.this_gene_model,
                                voucher_code
                            )
                        if self.aminoacids is not True:
                            sequence = self.degenerate(sequence, ThisGeneAndPartition.this_gene_model)

                        gene_codes_and_lengths[ThisGeneAndPartition.this_gene] = len(sequence)

                        tmp_seq = sequence.replace('-', '')
                        if len(tmp_seq) < 1:
                            out = ['\n[{}]'.format(voucher_code)]
                        else:
                            if self.outgroup != '' and self.outgroup not in voucher_code:
                                out += [line[0].ljust(55, ' ') + sequence]
                            elif self.outgroup == '':
                                out += [line[0].ljust(55, ' ') + sequence]

        out += ['\n;\nproc/;']

        number_chars = 0
        for k in gene_codes_and_lengths.keys():
            number_chars += gene_codes_and_lengths[k]

        header = ['nstates dna;\nxread']
        header += [str(number_chars) + ' ' + str(self.number_taxa - len(self.vouchers_to_drop))]

        out = header + out
        dataset_str = '\n'.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str
