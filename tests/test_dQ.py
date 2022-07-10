import unittest
from dQ.operation import Diann

class Test_dQ(unittest.TestCase):
    def test_dq_class(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_single_dq_class(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\DIA-NN_DIA-PASEF", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\DIA-NN_DIA-PASEF\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_single_microglia(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_LT-Micorglia", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_LT-Micorglia\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_CBQCA(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_CBQCA Test", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_CBQCA Test\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_CBQCA_gaf(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\20220419 - SCP IDs", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\20220419 - SCP IDs\temp")

    def test_stuff(self):
        d = Diann(r"C:\Users\Toan Phung\Desktop\deploy\data\dq\17f7d003-2a5e-4940-9016-e62c287b7747\data", r"C:\Users\Toan Phung\Desktop\deploy\data\dq\17f7d003-2a5e-4940-9016-e62c287b7747\data2")
if __name__ == '__main__':
    unittest.main()
