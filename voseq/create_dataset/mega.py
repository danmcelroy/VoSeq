from collections import namedtuple
from collections import OrderedDict

from .dataset import Dataset
from public_interface.models import Genes


class CreateMEGA(Dataset):
    def convert_lists_to_dataset(self, partitions):
        if self.partition_by_positions != 'ONE':
            self.errors += ['Cannot produce MEGA dataset with codon positions'
                            ' in different partitions.']
            return ''

        sequences_dict = OrderedDict()
        gene_code = ''
        for partition in partitions:
            for line in partition:
                parsed_line = self.parse_line(line)
                if parsed_line.gene_code != '':
                    gene_code = parsed_line.gene_code

                if parsed_line.taxon:
                    taxon = parsed_line.taxon.replace('?', '')
                    voucher_code = taxon.split('_')[0]
                    sequence = parsed_line.sequence

                    if self.partition_by_positions != 'ONE' and self.translations is True:
                        self.warnings.append('Cannot degenerate codons if they go to different partitions.')
                        continue

                    if self.aminoacids is True:
                        sequence = self.translate_this_sequence(
                            sequence,
                            gene_code,
                            voucher_code,
                        )

                    if self.aminoacids is not True and gene_code != '' and \
                            self.translations is True:
                        sequence = self.degenerate(sequence, gene_code)

                    if taxon not in sequences_dict:
                        sequences_dict[taxon] = ''
                    sequences_dict[taxon] += sequence

        out = ''
        for taxon, sequences in sequences_dict.items():
            out += '#{}\n{}\n'.format(taxon, sequences)

        dataset_str = '#MEGA\n!TITLE title;\n\n' + out.strip()
        self.save_dataset_to_file(dataset_str)
        return dataset_str

    def parse_line(self, line):
        """
        Takes an item from our list of partitions and tries to find the
        gene_code used. Also parses our item and returns items to be used for
        constructing our line to be used in the dataset.

        :param line:
        :param this_gene:
        :return: named tuple: tuple.gene_code, tuple.taxon, tuple.sequence
        """
        ParsedLine = namedtuple('ParsedLine', ['gene_code', 'taxon', 'sequence'])

        if line.startswith('\n['):
            this_gene_partition = line.replace('[', '').replace(']', '').strip()
            this_gene = Genes.objects.filter(
                gene_code=this_gene_partition).values()
            try:
                parsed_line = ParsedLine(this_gene[0], None, None)
            except KeyError:
                parsed_line = ParsedLine('', None, None)
            return parsed_line

        else:
            split_line = line.split(' ')
            voucher_code = split_line[0]
            sequence = split_line[-1]
            parsed_line = ParsedLine('', voucher_code, sequence)
            return parsed_line
