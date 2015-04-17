from collections import namedtuple
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import flatten_taxon_names_dict
from core.utils import translate_to_protein
from core.utils import strip_question_marks
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers
from .dataset import Dataset


class CreateFasta(Dataset):
    pass


class CreateGenbankFasta(Dataset):
    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
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
                                aa_sequence, warning = translate_to_protein(this_gene_model, sequence, '', voucher_code, self.file_format)
                                if warning != '':
                                    self.warnings.append(warning)

                            if aa_sequence.strip() == '':
                                self.warnings.append("Sequence for %s %s was empty" % (voucher_code, this_gene))

                            sequence = strip_question_marks(sequence)[0]
                            out += [line[0] + '\n' + sequence.replace('?', 'N') + '\n']
                            aa_out += [line[0] + '\n' + aa_sequence + '\n']

        dataset_str = ''.join(out)
        aa_dataset_str = ''.join(aa_out)
        self.save_dataset_to_file(dataset_str)
        self.save_aa_dataset_to_file(aa_dataset_str)
        return dataset_str, aa_dataset_str


class CreatePhylip(Dataset):
    def get_charset_block(self, gene_codes_and_lengths):
        charset_block = []

        bp_count_start = 0
        bp_count_end = 0
        self.gene_codes.sort()
        for gene in gene_codes_and_lengths:
            bp_count_end += gene_codes_and_lengths[gene]
            line = 'DNA, ' + gene + ' = ' + str(
                bp_count_start + 1) + '-' + str(bp_count_end)
            bp_count_start += gene_codes_and_lengths[gene]
            charset_block.append(line)
        self.charset_block = "\n".join(charset_block)

    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
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

                        gene_codes_and_lengths[ThisGeneAndPartition.this_gene] = len(sequence)

                        if partitions_incorporated == 1:
                            out += [line[0].ljust(55, ' ') + sequence + '\n']
                        else:
                            out += [' ' * 55 + sequence + '\n']

        number_chars = 0
        for k, v in gene_codes_and_lengths.items():
            number_chars += gene_codes_and_lengths[k]

        header = [
            str(self.number_taxa - len(self.vouchers_to_drop)) + ' ' + str(number_chars),
        ]

        out = header + out

        self.get_charset_block(gene_codes_and_lengths)
        dataset_str = ''.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str


class CreateTNT(Dataset):
    def convert_lists_to_dataset(self, partitions):
        """
        Overriden method from base clase in order to add headers and footers depending
        on needed dataset.
        """
        self.get_number_of_genes_for_taxa(partitions)
        self.get_number_chars_from_partition_list(partitions)

        out = 'nstates dna;\nxread\n'
        out += str(self.number_chars) + ' ' + str(self.number_taxa - len(self.vouchers_to_drop))

        outgroup_sequences = []
        for partition in partitions:
            for i in partition:
                voucher = i.split(' ')[0]
                if voucher not in self.vouchers_to_drop:
                    if self.outgroup in voucher:
                        outgroup_sequences.append(i)
                        continue

        partition_count = 0
        for partition in partitions:
            for i in partition:
                voucher = i.split(' ')[0]
                if voucher == '\n[&dna]':
                    out += '\n' + i
                    if self.outgroup != '':
                        out += '\n' + outgroup_sequences[partition_count]
                        partition_count += 1
                elif voucher not in self.vouchers_to_drop:
                    if self.outgroup != '' and self.outgroup not in voucher:
                        out += '\n' + i
                    elif self.outgroup == '':
                        out += '\n' + i

        out += '\n;\nproc/;'
        dataset_str = out.strip()
        self.save_dataset_to_file(dataset_str)
        return dataset_str


class CreateNEXUS(Dataset):
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

        gene_models = Genes.objects.all().values()
        gene_codes_and_lengths = dict()

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

                    gene_codes_and_lengths[ThisGeneAndPartition.this_gene] = len(sequence)

                    out += [line[0].ljust(55, ' ') + sequence]

        number_chars = 0
        for k, v in gene_codes_and_lengths.items():
            number_chars += gene_codes_and_lengths[k]

        header = [
            '#NEXUS\n',
            'BEGIN DATA;',
            'DIMENSIONS NTAX=' + str(self.number_taxa - len(self.vouchers_to_drop)) + ' NCHAR=' + str(number_chars) + ';',
            'FORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;',
            'MATRIX',
        ]

        out = header + out

        out += [';\nEND;']
        out += ['\nbegin mrbayes;']
        out += self.get_charset_block()
        out += self.get_partitions_block()
        out += self.get_final_block()
        dataset_str = '\n'.join(out)
        self.save_dataset_to_file(dataset_str)
        return dataset_str


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
        self.aminoacids = cleaned_data['aminoacids']
        self.codon_positions = cleaned_data['positions']
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
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

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
                                self.file_format)
            fasta_dataset = fasta.from_seq_objs_to_dataset()
            self.warnings += fasta.warnings
            self.dataset_file = fasta.dataset_file
            return fasta_dataset

        if self.file_format == 'PHY':
            phy = CreatePhylip(self.codon_positions, self.partition_by_positions,
                               self.seq_objs, self.gene_codes, self.voucher_codes,
                               self.file_format, self.outgroup, self.voucher_codes_metadata,
                               self.minimum_number_of_genes, self.aminoacids)
            phylip_dataset = phy.from_seq_objs_to_dataset()
            self.warnings += phy.warnings
            self.dataset_file = phy.dataset_file
            self.charset_block = phy.charset_block
            return phylip_dataset

        if self.file_format == 'TNT':
            tnt = CreateTNT(self.codon_positions, self.partition_by_positions,
                            self.seq_objs, self.gene_codes, self.voucher_codes,
                            self.file_format, self.outgroup, self.voucher_codes_metadata,
                            self.minimum_number_of_genes)
            tnt_dataset = tnt.from_seq_objs_to_dataset()
            self.warnings += tnt.warnings
            self.dataset_file = tnt.dataset_file
            return tnt_dataset

        if self.file_format == 'NEXUS':
            nexus = CreateNEXUS(self.codon_positions, self.partition_by_positions,
                                self.seq_objs, self.gene_codes, self.voucher_codes,
                                self.file_format, self.outgroup, self.voucher_codes_metadata,
                                self.minimum_number_of_genes, self.aminoacids)
            nexus_dataset = nexus.from_seq_objs_to_dataset()
            self.warnings += nexus.warnings
            self.dataset_file = nexus.dataset_file
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
        self.warnings += ['Could not find sequences for voucher %s' % i for i in vouchers_not_found]
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
