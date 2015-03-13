from django.test import TestCase
from django.core.management import call_command

from core import utils
from public_interface.models import Genes
from public_interface.models import Sequences


class TestCore(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

    def test_strip_question_marks_n(self):
        seq = '?NNNNNACTACGATGCRGCAST'
        result = utils.strip_question_marks(seq)
        expected = ('ACTACGATGCRGCAST', 6)
        self.assertEqual(expected, result)

    def test_translation_to_protein(self):
        """Catch exceptions when input has invalid codons due to ?"""
        gene_model = Genes.objects.filter(gene_code='COI').values()[0]
        sequence_model = Sequences.objects.get(gene_code='COI', code='CP100-10')
        seq_description = 'seq_description'
        seq_id = 'seq_id'
        expected = """
>seq_id seq_description
WAGMIGTSLSLIIRTELGNPSFLIGDDQIYNTIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNYMSFWLLPPSLILLISSSIVENGAGTGWTVYPPLSSNIAHSGASVDLAIFSLHLAGISSILGAINFITTIINMRINNMSYDQMPLFVWAVGITALLLLLSLPVLAGAITMLLTDRNLNTSFFDSCGGGD
"""
        results = utils.translate_to_protein(gene_model, sequence_model, seq_description, seq_id)
        self.assertEqual(expected.lstrip(), results)

    def test_translation_to_protein_invalid_codons(self):
        """Catch exceptions when input has invalid codons due to ?"""
        gene_model = Genes.objects.filter(gene_code='COI').values()[0]
        sequence_model = Sequences.objects.get(gene_code='COI', code='CP100-10')
        seq = sequence_model.sequences
        new_seq = []
        count = 1
        for i in seq:
            if count % 10 == 0:
                i = '?'
            new_seq.append(i)
            count += 1
        sequence_model.sequences = ''.join(new_seq)
        seq_description = 'seq_description'
        seq_id = 'seq_id'
        expected = "Error Codon 'G?C' is invalid"
        results = utils.translate_to_protein(gene_model, sequence_model, seq_description, seq_id)
        self.assertEqual(expected, results)
