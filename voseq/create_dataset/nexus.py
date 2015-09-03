import collections

from .dataset import Dataset
from public_interface.models import Genes


class CreateNEXUS(Dataset):
    def get_partitions_block(self):
        line = 'partition GENES = '
        number_of_codons = self.get_number_of_codons()
        if self.partition_by_positions == 'EACH':
            line += str(len(self.gene_codes_and_lengths) * number_of_codons)
        elif self.partition_by_positions == '1st2nd_3rd':
            line += str(len(self.gene_codes_and_lengths) * number_of_codons)
        else:
            line += str(len(self.gene_codes_and_lengths))
        line += self.make_partition_line()
        line += '\nset partition = GENES;\n'
        return [line]

    def get_number_of_codons(self):
        if 'ALL' in self.codon_positions and self.partition_by_positions == 'ONE':
            return 1
        elif 'ALL' in self.codon_positions and self.partition_by_positions == '1st2nd_3rd':
            return 2
        elif 'ALL' in self.codon_positions and self.partition_by_positions == 'EACH':
            return 3
        elif '1st' in self.codon_positions and '2nd' in self.codon_positions and self.partition_by_positions == 'ONE':
            return 1
        elif '1st' in self.codon_positions and '2nd' in self.codon_positions and self.partition_by_positions == '1st2nd_3rd':
            return 1
        elif '1st' in self.codon_positions and '2nd' in self.codon_positions and self.partition_by_positions == 'EACH':
            return 2
        else:
            return 1

    def make_partition_line(self):
        if self.partition_by_positions == 'ONE':
            if 'ALL' not in self.codon_positions and len(self.codon_positions) == 1:
                return self.build_gene_line_for_one_codon_position()
            elif len(self.codon_positions) == 2:
                return self.build_gene_line_for_two_codon_positions()
            else:
                return ': ' + ', '.join([i for i in self.gene_codes_and_lengths]) + ';\n'

        if self.partition_by_positions == 'EACH':
            if 'ALL' not in self.codon_positions and len(self.codon_positions) == 1:
                return self.build_gene_line_for_one_codon_position()
            elif len(self.codon_positions) == 2:
                return self.build_gene_line_for_two_codon_positions_partitioned_each()
            else:
                out = []
                for i in self.gene_codes_and_lengths:
                    out += ['{}_pos1'.format(i), '{}_pos2'.format(i), '{}_pos3'.format(i)]
                return ': ' + ', '.join(out) + ';\n'
        if self.partition_by_positions == '1st2nd_3rd':
            if 'ALL' not in self.codon_positions and len(self.codon_positions) == 1:
                return self.build_gene_line_for_one_codon_position()
            elif len(self.codon_positions) == 2:
                return self.build_gene_line_for_two_codon_positions()
            else:
                out = []
                for i in self.gene_codes_and_lengths:
                    out += ['{}_pos12'.format(i), '{}_pos3'.format(i)]
                return ': ' + ', '.join(out) + ';\n'

    def build_gene_line_for_one_codon_position(self):
        out = []
        for i in self.gene_codes_and_lengths:
            if '1st' in self.codon_positions:
                out += ['{}_pos1'.format(i)]
            elif '2nd' in self.codon_positions:
                out += ['{}_pos2'.format(i)]
            elif '3rd' in self.codon_positions:
                out += ['{}_pos3'.format(i)]
        return ': ' + ', '.join(out) + ';\n'

    def build_gene_line_for_two_codon_positions(self):
        out = []
        for i in self.gene_codes_and_lengths:
            if '1st' in self.codon_positions and '2nd' in self.codon_positions:
                out += ['{}_pos12'.format(i)]
        return ': ' + ', '.join(out) + ';\n'

    def build_gene_line_for_two_codon_positions_partitioned_each(self):
        out = []
        for i in self.gene_codes_and_lengths:
            if '1st' in self.codon_positions and '2nd' in self.codon_positions:
                out += ['{}_pos1'.format(i)]
                out += ['{}_pos2'.format(i)]
        return ': ' + ', '.join(out) + ';\n'

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
        Overridden method from base class in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        gene_models = Genes.objects.all().values()
        gene_codes_and_lengths = collections.OrderedDict()

        partitions_incorporated = 0

        out = []

        for partition in partitions:
            for i in partition:
                voucher_code = i.split(' ')[0]
                if voucher_code.startswith('\n'):
                    ThisGeneAndPartition = self.get_gene_for_current_partition(
                        gene_models, out, partitions_incorporated, voucher_code
                    )
                    out = ThisGeneAndPartition.out
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
                        if tmp_seq:
                            out += [line[0].ljust(55, ' ') + sequence]
                        else:
                            out = ['\n[{}]'.format(voucher_code)]

        number_chars = 0
        for k in gene_codes_and_lengths.keys():
            number_chars += gene_codes_and_lengths[k]

        header = [
            '#NEXUS\n',
            'BEGIN DATA;',
            'DIMENSIONS NTAX=' + str(self.number_taxa - len(self.vouchers_to_drop)) + ' NCHAR=' + str(number_chars) + ';',
            'FORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;',
            'MATRIX',
        ]

        if self.aminoacids is True:
            self.gene_codes_and_lengths = gene_codes_and_lengths

        self.make_charset_block(self.gene_codes_and_lengths)

        out = header + out

        out += [';\nEND;']
        out += ['\nbegin mrbayes;']
        out += [self.charset_block]
        out += self.get_partitions_block()
        out += self.get_final_block()
        dataset_str = '\n'.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str
