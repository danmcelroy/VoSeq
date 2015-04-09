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
        self.warnings = []
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

    def build_fasta_seq_components(self, code, gene, sequence_model,
                                   voucher_models):
        for v in voucher_models:
            if v['code'] == code:
                seq_id = v['genus'] + '_' + v['species'] + '_' + code
                seq_description = '[org=' + v['genus'] + ' ' + v[
                    'species'] + ']'
                seq_description += ' [Specimen-voucher=' + code + ']'
                seq_description += ' [note=' + gene[
                    'description'] + ' gene, partial cds.]'
                seq_description += ' [Lineage=]'
                seq_seq = sequence_model.sequences

        return seq_description, seq_id, seq_seq

    def get_datasets(self):
        """Queries sequences and creates FASTA, protein strings and list of
        items with accession number (code, gene_code, accession).
        """
        sequence_models = Sequences.objects.all()
        voucher_models = Vouchers.objects.all().values('genus', 'species', 'code')
        gene_models = Genes.objects.all().values()

        vouchers = self.get_vouchers_from_voucher_models(voucher_models)
        genes = self.get_genes_from_gene_models(gene_models)

        for sequence_model in sequence_models:
            code = sequence_model.code_id
            gene_code = sequence_model.gene_code
            gene = genes[gene_code.lower()]
            if code in self.voucher_codes and gene_code in self.gene_codes:
                if sequence_model.accession.strip() != '':
                    self.items_with_accession.append(
                        {
                            'voucher_code': code,
                            'gene_code': gene_code,
                            'accession': sequence_model.accession,
                        },
                    )
                else:
                    seq_description, seq_id, seq_seq = self.build_fasta_seq_components(
                        code, gene, sequence_model, voucher_models)

                    # # DNA sequences
                    seq_seq = utils.strip_question_marks(seq_seq)[0]
                    if '?' in seq_seq or 'N' in seq_seq.upper():
                        seq_obj = Seq(seq_seq, IUPAC.ambiguous_dna)
                    else:
                        seq_obj = Seq(seq_seq, IUPAC.unambiguous_dna)

                    self.fasta += '>' + seq_id + ' ' + seq_description + '\n'
                    self.fasta += str(seq_obj) + '\n'

                    protein, warning = utils.translate_to_protein(
                        gene,
                        sequence_model.sequences,
                        seq_description,
                        seq_id,
                        'FASTA',
                    )
                    if warning != '':
                        self.warnings.append(warning)

                    if not protein.startswith('Error'):
                        self.protein += protein
                    else:
                        self.warnings.append("Could not translate %s: %s" % (seq_id, protein))

        with open(self.fasta_file, 'w') as handle:
            handle.write(self.fasta)

        with open(self.protein_file, 'w') as handle:
            handle.write(self.protein)

    def get_vouchers_from_voucher_models(self, voucher_models):
        vouchers = dict()
        for voucher_model in voucher_models:
            code = voucher_model['code']
            vouchers[code] = voucher_model
        return vouchers

    def get_genes_from_gene_models(self, gene_models):
        genes = dict()
        for gene_model in gene_models:
            gene_code = gene_model['gene_code'].lower()
            genes[gene_code] = gene_model
        return genes

    def make_guid(self):
        return uuid.uuid4().hex
