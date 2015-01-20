from django.test import TestCase
from django.core.management import call_command

from genbank_fasta.utils import Results


class ResultTest(TestCase):
    def setUp(self):
        args = []
        opts = {'dumpfile': 'test_db_dump.xml', 'verbosity': 0}
        cmd = 'migrate_db'
        call_command(cmd, *args, **opts)

        self.maxDiff = None

    def test_fasta(self):
        res = Results(['CP100-10', 'CP100-11'], ['COI'])
        res.get_datasets()

        expected = '>Melitaea_diamina_CP100-10 [org=Melitaea diamina] [Specimen-voucher=CP100-10] [note=cytochrome oxidase c subunit I gene, partial cds.] [Lineage=]\nTGAGCCGGTATAATTGGTACATCCCTAAGTCTTATTATTCGAACCGAATTAGGAAATCCTAGTTTTTTAATTGGAGATGATCAAATTTATAATACCATTGTAACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATGCCAATTATAATTGGAGGATTTGGTAATTGACTTGTACCATTAATATTGGGAGCCCCAGATATAGCTTTCCCCCGAATAAATTATATAAGATTTTGATTATTGCCTCCATCCTTAATTCTTTTAATTTCAAGTAGAATTGTAGAAAATGGGGCAGGAACTGGATGAACAGTTTACCCCCCACTTTCATCTAATATTGCCCATAGAGGAGCTTCAGTGGATTTAGCTATTTTTTCTTTACATTTAGCTGGGATTTCCTCTATCTTAGGAGCTATTAATTTTATTACTACAATTATTAATATACGAATTAATAATATATCTTATGATCAAATACCTTTATTTGTATGAGCAGTAGGAATTACAGCATTACTTCTCTTATTATCTTTACCAGTTTTAGCTGGAGCTATTACTATACTTTTAACGGATCGAAATCTTAATACCTCATTTTTTGATTCCTGCGGAGGAGGAGATCC\n'
        self.assertEqual(expected, res.fasta)

    def test_protein(self):
        res = Results(['CP100-10', 'CP100-11'], ['COI'])
        res.get_datasets()

        expected = '>Melitaea_diamina_CP100-10 [org=Melitaea diamina] [Specimen-voucher=CP100-10] [note=cytochrome oxidase c subunit I gene, partial cds.] [Lineage=]\nWAGMIGTSLSLIIRTELGNPSFLI\n'
        self.assertEqual(expected, res.protein)
