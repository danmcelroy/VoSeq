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
        results = utils.translate_to_protein(gene_model, sequence_model.sequences, seq_description, seq_id)
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
        expected = ''
        results = utils.translate_to_protein(gene_model, sequence_model.sequences, seq_description, seq_id)
        self.assertEqual(expected, results)

    def test_flatten_taxon_names_dict(self):
        dictionary = {'code': 'CP100-10', 'orden': 'Lepidoptera',
                      'genus': 'Danaus', 'species': 'sp. 5'}
        expected = 'CP100-10_Lepidoptera_Danaus_sp._5'
        results = utils.flatten_taxon_names_dict(dictionary)
        self.assertEqual(expected, results)

    def test_gapped_translation(self):
        sequence = 'ATG---GCCATTGTAATGGGCCGG'
        expected = 'M?AIVMGR'
        result = utils.gapped_translation(sequence)
        self.assertEqual(expected, result)

    def test_getting_gap_indexes1(self):
        sequence = 'ATG---GCCATTGTAATGGGCCGG'
        expected = [1]
        result = utils.get_gap_indexes(sequence)
        self.assertEqual(expected, result)

    def test_getting_gap_indexes2(self):
        sequence = 'ATG---GCCATT---GTAATGGGCCGG'
        expected = [1, 4]
        result = utils.get_gap_indexes(sequence)
        self.assertEqual(expected, result)

    def test_getting_gap_indexes3(self):
        """Dont mind if the are invalid codons. They will rejected by the
        translation tool of Biopython.
        """
        sequence = 'ATG---GCC-TTGTAATGGGCCGG'
        expected = [1]
        result = utils.get_gap_indexes(sequence)
        self.assertEqual(expected, result)
