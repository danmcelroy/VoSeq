from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from core.utils import flatten_taxon_names_dict
from core.utils import chain_and_flatten
from public_interface.models import Genes
from public_interface.models import Sequences
from public_interface.models import Vouchers


class CreateFasta(object):
    def __init__(self, codon_positions, partition_by_positions, seq_objs, gene_codes):
        self.codon_positions = codon_positions
        self.partition_by_positions = partition_by_positions
        self.seq_objs = seq_objs
        self.gene_codes = gene_codes
        self.reading_frames = self.get_reading_frames()
        self.warnings = []

    def get_reading_frames(self):
        """

        :return: dict of gene_code: reading_frame. If not found, flag warning.
        """
        reading_frames = dict()
        genes = Genes.objects.all().values('gene_code', 'reading_frame')
        for gene in genes:
            gene_code = gene['gene_code'].lower()
            if gene_code in self.gene_codes:
                reading_frames[gene_code] = gene['reading_frame']
        return reading_frames

    def split_sequence_in_codon_positions(self, gene_code, seq):
        """Puts the sequence in frame, by deleting base pairs at the beginning
        of the sequence if the reading frame is not 1.

        Retuns tuple of nucleotides based on codon positions.

        :param gene_code: as lower case
        :param seq: as BioPython seq object.
        :return: tuple of seq strings.

        Example:
            If reading frame is 2: ATGGGG becomes TGGGG. Then the sequence is
            processed to extract the codon positions requested by the user.

        """
        reading_frame = int(self.reading_frames[gene_code.lower()]) - 1
        seq = seq[reading_frame:]

        # This is the BioPython way to get codon positions
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc19
        first_position = seq[0::3]
        second_position = seq[1::3]
        third_position = seq[2::3]
        return first_position, second_position, third_position

    def convert_lists_to_dataset(self, partitions):
        out = ''
        for i in partitions:
            out += '\n'
            out += '\n'.join(i)
        return out.strip()

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
                        seq_str = '>' + this_gene + codon_descriptions[i] + '\n--------------------'
                        partition_list[i].append(seq_str)

                codons = self.split_sequence_in_codon_positions(this_gene,
                                                                seq_record.seq)
                for i in range(len(codon_pos)):
                    seq_str = '>' + seq_record.id + '\n' + str(codons[codon_pos[i]])
                    partition_list[i].append(seq_str)
        return partition_list

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

                    seq_str = '>' + this_gene + '\n' + '--------------------'
                    partition_list[0].append(seq_str)

                codons = self.split_sequence_in_codon_positions(this_gene,
                                                                seq_record.seq)

                seq_str = '>' + seq_record.id + '\n'
                seq_str += str(chain_and_flatten([codons[i] for i in codon_pos]))

                partition_list[0].append(seq_str)
        return partition_list

    def from_seq_objs_to_fasta(self):
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
            partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                'ALL' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_one_partition(['1st', '2nd'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '3rd' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['1st', '3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '3rd' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_one_partition(['1st', '3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '2nd' in self.codon_positions and '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['2nd', '3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '2nd' in self.codon_positions and '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                'ONE' in self.partition_by_positions:
            partition_list = self.get_codons_in_one_partition(['2nd', '3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if ('ALL' in self.codon_positions or
                ('1st' in self.codon_positions and '2nd' in self.codon_positions and '3rd' in self.codon_positions)) \
                and 'EACH' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['1st', '2nd', '3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if ('ALL' in self.codon_positions or
                ('1st' in self.codon_positions and '2nd' in self.codon_positions and '3rd' in self.codon_positions)) \
                and 'ONE' in self.partition_by_positions:
            partition_list = ([],)
            for gene_code in self.seq_objs:
                this_gene = None
                for seq_record in self.seq_objs[gene_code]:
                    if self.reading_frames[gene_code] is None:
                        self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                             "it cannot be included in your dataset." % gene_code)
                        continue

                    if this_gene is None:
                        this_gene = seq_record.name

                        seq_str = '>' + this_gene + '\n' + '--------------------'
                        partition_list[0].append(seq_str)

                    seq_str = '>' + seq_record.id + '\n' + str(seq_record.seq)
                    partition_list[0].append(seq_str)
            return self.convert_lists_to_dataset(partition_list)

        if 'ALL' in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = ([], [],)

            for gene_code in self.seq_objs:
                if self.reading_frames[gene_code] is None:
                    self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                         "it cannot be included in your dataset." % gene_code)
                    continue

                this_gene = None
                for seq_record in self.seq_objs[gene_code]:
                    if self.reading_frames[gene_code] is None:
                        self.warnings.append("Reading frame for gene %s hasn't been specified so "
                                             "it cannot be included in your dataset." % gene_code)
                        continue

                    if this_gene is None:
                        this_gene = seq_record.name

                        seq_str = '>' + this_gene + '_1st_2nd_codons\n' + '--------------------'
                        partition_list[0].append(seq_str)

                        seq_str = '>' + this_gene + '_3rd_codon\n' + '--------------------'
                        partition_list[1].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(chain_and_flatten([codons[0], codons[1]]))
                    partition_list[0].append(seq_str)

                    seq_str = '>' + seq_record.id + '\n' + str(codons[2])
                    partition_list[1].append(seq_str)
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['1st'])
            return self.convert_lists_to_dataset(partition_list)

        if '2nd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['2nd'])
            return self.convert_lists_to_dataset(partition_list)

        if '3rd' in self.codon_positions and \
                '1st' not in self.codon_positions and \
                '2nd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = self.get_codons_in_each_partition(['3rd'])
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' not in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = ([],)
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
                        partition_list[0].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(chain_and_flatten([codons[0], codons[1]]))
                    partition_list[0].append(seq_str)
            return self.convert_lists_to_dataset(partition_list)

        if '1st' in self.codon_positions and '2nd' in self.codon_positions and \
                '3rd' in self.codon_positions and \
                '1st2nd_3rd' in self.partition_by_positions:
            partition_list = ([], [],)
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
                        partition_list[0].append(seq_str)

                        seq_str = '>' + this_gene + '_3rd_codon\n' + '--------------------'
                        partition_list[1].append(seq_str)

                    codons = self.split_sequence_in_codon_positions(this_gene, seq_record.seq)

                    seq_str = '>' + seq_record.id + '\n' + str(chain_and_flatten([codons[0], codons[1]]))
                    partition_list[0].append(seq_str)

                    seq_str = '>' + seq_record.id + '\n' + str(codons[2])
                    partition_list[1].append(seq_str)
            return self.convert_lists_to_dataset(partition_list)


class CreateDataset(object):
    """
    Accept form input to create a dataset in several formats, codon positions,
    for list of codes and genes. Also takes into account the vouchers passed as
    taxonset.

    Attributes:
        ``dataset_str``: output dataset to pass to users.

    """
    def __init__(self, cleaned_data):
        print(">>>>>>_init", cleaned_data)
        self.errors = []
        self.seq_objs = dict()
        self.codon_positions = cleaned_data['positions']
        self.partition_by_positions = cleaned_data['partition_by_positions']
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.taxon_names = cleaned_data['taxon_names']
        self.warnings = []
        self.dataset_str = self.create_dataset()

    def create_dataset(self):
        self.voucher_codes = get_voucher_codes(self.cleaned_data)
        self.gene_codes = get_gene_codes(self.cleaned_data)
        self.create_seq_objs()

        fasta = CreateFasta(self.codon_positions, self.partition_by_positions, self.seq_objs, self.gene_codes)
        fasta_dataset = fasta.from_seq_objs_to_fasta()
        self.warnings += fasta.warnings
        return fasta_dataset

    def create_seq_objs(self):
        """Generate a list of sequence objects. Also takes into account the
        genes passed as geneset.

        Returns:
            list of sequence objects as produced by BioPython.

        """
        our_taxon_names = self.get_taxon_names_for_taxa()

        all_seqs = Sequences.objects.all().values('code_id', 'gene_code', 'sequences').order_by('code_id')
        for s in all_seqs:
            code = s['code_id'].lower()
            gene_code = s['gene_code'].lower()
            if code in self.voucher_codes and gene_code in self.gene_codes:
                seq = Seq(s['sequences'])
                seq_obj = SeqRecord(seq)
                seq_obj.id = flatten_taxon_names_dict(our_taxon_names[code])
                if 'GENECODE' in self.taxon_names:
                    seq_obj.id += '_' + gene_code
                seq_obj.name = gene_code

                if gene_code not in self.seq_objs:
                    self.seq_objs[gene_code] = []
                self.seq_objs[gene_code].append(seq_obj)

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
            code = voucher['code'].lower()
            if code in self.voucher_codes:
                obj = dict()
                for taxon_name in self.taxon_names:
                    if taxon_name != 'GENECODE':
                        taxon_name = taxon_name.lower()
                        obj[taxon_name] = voucher[taxon_name]
                vouchers_with_taxon_names[code] = obj

        return vouchers_with_taxon_names
