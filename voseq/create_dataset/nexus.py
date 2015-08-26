import collections

from .dataset import Dataset
from public_interface.models import Genes


class CreateNEXUS(Dataset):
    def get_partitions_block(self):
        line = 'partition GENES = ' + str(len(self.gene_codes_and_lengths))
        line += self.make_partition_line()
        line += '\nset partition = GENES;\n'
        return [line]

    def make_partition_line(self):
        if self.partition_by_positions == 'ONE':
            return ': ' + ', '.join([i for i in self.gene_codes_and_lengths]) + ';\n'
        if self.partition_by_positions == 'EACH':
            out = []
            for i in self.gene_codes_and_lengths:
                out += ['{}_pos1'.format(i), '{}_pos2'.format(i), '{}_pos3'.format(i)]
            return ': ' + ', '.join(out) + ';\n'
        if self.partition_by_positions == '1st2nd_3rd':
            out = []
            for i in self.gene_codes_and_lengths:
                out += ['{}_pos12'.format(i), '{}_pos3'.format(i)]
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
                        if len(tmp_seq) < 1:
                            out = ['\n[{}]'.format(voucher_code)]
                        else:
                            out += [line[0].ljust(55, ' ') + sequence]

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
