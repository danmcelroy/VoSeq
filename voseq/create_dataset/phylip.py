from .dataset import Dataset
from public_interface.models import Genes


class CreatePhylip(Dataset):
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
        if self.partition_by_positions == 'ONE':
            if 'ALL' in self.codon_positions:
                line = 'DNA, {} = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = 'DNA, {}_pos1 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = 'DNA, {}_pos2 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = 'DNA, {}_pos3 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 2 and \
                    '1st' in self.codon_positions and \
                    '2nd' in self.codon_positions:
                line = 'DNA, {}_pos12 = {}-{}'.format(gene, count_start + 1, count_end)
                return line

        elif self.partition_by_positions == '1st2nd_3rd':
            if 'ALL' in self.codon_positions:
                line = ''
                if self.reading_frames[gene] == 1:
                    line = 'DNA, {}_pos12 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 1, count_end)
                    line += ', {}-{}\\3'.format(count_start + 2, count_end)
                    line += '\nDNA, {}_pos3 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 3, count_end)
                elif self.reading_frames[gene] == 2:
                    line = 'DNA, {}_pos12 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 2, count_end)
                    line += ', {}-{}\\3'.format(count_start + 3, count_end)
                    line += '\nDNA, {}_pos3 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 1, count_end)
                elif self.reading_frames[gene] == 3:
                    line = 'DNA, {}_pos12 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 3, count_end)
                    line += ', {}-{}\\3'.format(count_start + 1, count_end)
                    line += '\nDNA, {}_pos3 = '.format(gene)
                    line += '{}-{}\\3'.format(count_start + 2, count_end)
                return line

            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = 'DNA, {}_pos1 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = 'DNA, {}_pos2 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = 'DNA, {}_pos3 = {}-{}'.format(gene, count_start + 1, count_end)
                return line

        elif self.partition_by_positions == 'EACH':
            if 'ALL' in self.codon_positions:
                line = ''
                if self.reading_frames[gene] == 1:
                    line = 'DNA, {}_pos1 = {}-{}\\3'.format(gene, count_start + 1, count_end)
                    line += '\nDNA, {}_pos2 = {}-{}\\3'.format(gene, count_start + 2, count_end)
                    line += '\nDNA, {}_pos3 = {}-{}\\3'.format(gene, count_start + 3, count_end)
                elif self.reading_frames[gene] == 2:
                    line = 'DNA, {}_pos1 = {}-{}\\3'.format(gene, count_start + 2, count_end)
                    line += '\nDNA, {}_pos2 = {}-{}\\3'.format(gene, count_start + 3, count_end)
                    line += '\nDNA, {}_pos3 = {}-{}\\3'.format(gene, count_start + 1, count_end)
                elif self.reading_frames[gene] == 3:
                    line = 'DNA, {}_pos1 = {}-{}\\3'.format(gene, count_start + 3, count_end)
                    line += '\nDNA, {}_pos2 = {}-{}\\3'.format(gene, count_start + 1, count_end)
                    line += '\nDNA, {}_pos3 = {}-{}\\3'.format(gene, count_start + 2, count_end)
                return line

            if len(self.codon_positions) == 1 and '1st' in self.codon_positions:
                line = 'DNA, {}_pos1 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '2nd' in self.codon_positions:
                line = 'DNA, {}_pos2 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 1 and '3rd' in self.codon_positions:
                line = 'DNA, {}_pos3 = {}-{}'.format(gene, count_start + 1, count_end)
                return line
            if len(self.codon_positions) == 2 and \
                    '1st' in self.codon_positions and \
                    '2nd' in self.codon_positions:
                line = 'DNA, {}_pos12 = {}-{}'.format(gene, count_start + 1, count_end)
                return line

    def convert_lists_to_dataset(self, partitions):
        """
        Overridden method from base class in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)
        out = []

        gene_models = Genes.objects.all().values()

        gene_codes_and_lengths = dict()

        partitions_incorporated = 0
        for partition in partitions:
            for i in partition:
                voucher_code = i.split(' ')[0]
                if voucher_code.startswith('\n'):
                    ThisGeneAndPartition = self.get_gene_for_current_partition(
                        gene_models, out, partitions_incorporated, voucher_code
                    )
                    partitions_incorporated = ThisGeneAndPartition.partitions_incorporated
                    out = ThisGeneAndPartition.out
                elif voucher_code not in self.vouchers_to_drop:
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

                        if partitions_incorporated == 1:
                            out += [line[0].ljust(55, ' ') + sequence + '\n']
                        else:
                            out += [' ' * 55 + sequence + '\n']

        number_chars = 0
        for k in gene_codes_and_lengths.keys():
            number_chars += gene_codes_and_lengths[k]

        header = [
            str(self.number_taxa - len(self.vouchers_to_drop)) + ' ' + str(number_chars),
        ]

        out = header + out

        self.make_charset_block(gene_codes_and_lengths)
        dataset_str = ''.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str
