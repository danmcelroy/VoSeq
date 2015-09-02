from .dataset import Dataset
from public_interface.models import Genes


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
