import unittest
from dQ.operation import Diann

class Test_dQ(unittest.TestCase):
    def test_dq_class(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_single_dq_class(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\DIA-NN_DIA-PASEF", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\DIA-NN_DIA-PASEF\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")
    def test_single_microglia(self):
        d = Diann(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_LT-Micorglia", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For_LT-Micorglia\temp", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\Mouse_UniprotSP_Cano+Iso_052021.fasta")

if __name__ == '__main__':
    unittest.main()
