import os

from public_interface.models import Genes
from create_dataset import dataset


class Results(dataset.Dataset):
    def __init__(self, *args, **kwargs):
        super(Results, self).__init__(*args, **kwargs)
        self.gene_codes_and_lengths = None
        self.number_taxa = len(self.voucher_codes)
        self.number_chars = None
        self.vouchers_to_drop = None
        self.charset_block = None
        self.fasta_file = os.path.join(self.cwd,
                                       'fasta_files',
                                       'fasta_' + self.guid + '.fasta',
                                       )
        self.protein_file = os.path.join(self.cwd,
                                         'fasta_files',
                                         'prot_' + self.guid + '.fasta',
                                         )

    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        out = [
            str(self.number_taxa - len(self.vouchers_to_drop)) + ' ' + str(self.number_chars),
        ]

        gene_models = Genes.objects.all().values()

        partitions_incorporated = 0
        for partition in partitions:
            for i in partition:
                voucher_code = i.split(' ')[0]
                if voucher_code.startswith('\n'):
                    this_gene = voucher_code.replace('[', '').replace(']', '').strip()
                    this_gene_model = get_gene_model_from_gene_id(this_gene, gene_models)
                    partitions_incorporated += 1
                    out += ['\n']
                elif voucher_code not in self.vouchers_to_drop:
                    line = i.split(' ')
                    if len(line) > 1:
                        sequence = line[-1]

                        if self.aminoacids is True:
                            if this_gene_model['genetic_code'] is None or this_gene_model['reading_frame'] is None:
                                self.warnings.append("Cannot translate gene %s sequences into aminoacids."
                                                     " You need to define reading_frame and/or genetic_code." % this_gene_model['gene_code'])
                            else:
                                sequence, warning = translate_to_protein(this_gene_model, sequence, '', voucher_code, self.file_format)
                                if warning != '':
                                    self.warnings.append(warning)

                        if partitions_incorporated == 1:
                            out += [line[0].ljust(55, ' ') + sequence + '\n']
                        else:
                            out += [' ' * 55 + sequence + '\n']

        self.get_charset_block()
        dataset_str = ''.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str
