from Bio.Blast.Applications import NcbiblastnCommandline

from django.shortcuts import render


def index(request, voucher_code, gene_code):
    return render(request, 'public_interface/base.html')


class BLAST(object):
    def __init__(self, type, voucher_code, gene_code):
        """
        Type of blast to do: local, full, remote

        :param type: local, full, remote.
        :param voucher_code:
        :param gene_code:
        """
        self.blast_type = type
        self.voucher_code = voucher_code
        self.gener_code = gene_code

    def do_blast(self):
        blastn_cline = NcbiblastnCommandline(query=self.query, db=self.db,
                                             evalue=0.001, outfmt=5, out="opuntia.xml")
        blastn_cline()
