from Bio.Blast import NCBIWWW

from blast_local.utils import BLAST


class BLASTNcbi(BLAST):
    """Handles duties related to blast against sequences in NCBI GenBank."""
    def do_blast(self):
        """Blasts against NCBI and saves returned XML file to local disk."""
        with open(self.query_file) as handle:
            fasta_string = handle.read()

        result_handle = NCBIWWW.qblast('blastn', 'nt', fasta_string)

        with open(self.output_file, 'w') as writer:
            writer.write(result_handle.read())
        result_handle.close()
