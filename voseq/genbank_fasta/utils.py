import os
import uuid

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from public_interface.models import Vouchers
from public_interface.models import Sequences
from public_interface.models import Genes
from core import utils


class Results(object):
    """Returns:

    Attributes:
        * List of items with accession numbers that will not be included in FASTA or protein files (code, gene_code, accession).
        * FASTA dataset as string.
        * Protein dataset as string.

    Usage:

        >>> res = Results(voucher_codes, gene_codes)
        >>> res.get_datasets()
        >>> fasta_dataset = res.fasta
        >>> protein_dataset = res.protein
        >>> items_with_accession = res.items_with_accession

    """
    def __init__(self, voucher_codes, gene_codes):
        self.voucher_codes = voucher_codes
        self.gene_codes = gene_codes
        self.items_with_accession = []
        self.guid = self.make_guid()
        self.fasta = ''
        self.protein = ''
        self.cwd = os.path.dirname(__file__)
        self.fasta_file = os.path.join(self.cwd,
                                       'fasta_files',
                                       'fasta_' + self.guid + '.fasta',
                                       )
        self.protein_file = os.path.join(self.cwd,
                                         'fasta_files',
                                         'prot_' + self.guid + '.fasta',
                                         )

    def translate_to_protein(self, gene_model, sequence_model, seq_description, seq_id):
        # # Protein sequences
        seq_seq, removed = utils.strip_question_marks(sequence_model.sequences)
        if int(gene_model.reading_frame) == 1:
            if removed % 3 == 0:
                start_translation = 0
            if removed % 3 == 1:
                start_translation = 2
            if removed % 3 == 2:
                start_translation = 1
        if int(gene_model.reading_frame) == 2:
            if removed % 3 == 0:
                start_translation = 1
            if removed % 3 == 1:
                start_translation = 0
            if removed % 3 == 2:
                start_translation = 2
        if int(gene_model.reading_frame) == 3:
            if removed % 3 == 0:
                start_translation = 2
            if removed % 3 == 1:
                start_translation = 1
            if removed % 3 == 2:
                start_translation = 0
        if '?' in seq_seq or 'N' in seq_seq.upper():
            seq_obj = Seq(seq_seq[start_translation:], IUPAC.ambiguous_dna)
        else:
            seq_obj = Seq(seq_seq[start_translation:], IUPAC.unambiguous_dna)
        prot_sequence = seq_obj.translate(table=gene_model.genetic_code)
        out = '>' + seq_id + ' ' + seq_description + '\n'
        out += str(prot_sequence) + '\n'
        return out

    def get_datasets(self):
        """Queries sequences and creates FASTA, protein strings and list of
        items with accession number (code, gene_code, accession).
        """
        for code in self.voucher_codes:
            try:
                v = Vouchers.objects.get(code=code)
            except Vouchers.DoesNotExist:
                continue

            for gene_code in self.gene_codes:
                if v:
                    try:
                        sequence_model = Sequences.objects.get(code=v, gene_code=gene_code)
                    except Sequences.DoesNotExist:
                        continue

                    try:
                        gene_model = Genes.objects.get(gene_code=gene_code)
                    except Genes.DoesNotExist:
                        continue

                    if sequence_model.accession.strip() != '':
                        self.items_with_accession.append(
                            {
                                'voucher_code': code,
                                'gene_code': gene_code,
                                'accession': sequence_model.accession,
                            },
                        )
                    else:
                        seq_id = v.genus + '_' + v.species + '_' + code
                        seq_description = '[org=' + v.genus + ' ' + v.species + ']'
                        seq_description += ' [Specimen-voucher=' + v.code + ']'
                        seq_description += ' [note=' + gene_model.description + ' gene, partial cds.]'
                        seq_description += ' [Lineage=]'
                        seq_seq = sequence_model.sequences

                        # # DNA sequences
                        seq_seq = utils.strip_question_marks(seq_seq)[0]
                        if '?' in seq_seq or 'N' in seq_seq.upper():
                            seq_obj = Seq(seq_seq, IUPAC.ambiguous_dna)
                        else:
                            seq_obj = Seq(seq_seq, IUPAC.unambiguous_dna)

                        self.fasta += '>' + seq_id + ' ' + seq_description + '\n'
                        self.fasta += str(seq_obj) + '\n'

                        self.protein += self.translate_to_protein(
                            gene_model,
                            sequence_model,
                            seq_description,
                            seq_id,
                        )

        with open(self.fasta_file, 'w') as handle:
            handle.write(self.fasta)

        with open(self.protein_file, 'w') as handle:
            handle.write(self.protein)

    def make_guid(self):
        return uuid.uuid4().hex
