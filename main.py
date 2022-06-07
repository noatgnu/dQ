import pandas as pd
import pandera as pa
from uniprot.parser import acc_regex

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    fasta_library = r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For Phospho Motif\Updated\Human_UniprotSP_Cano+Iso_052021.fasta"
    seqs = {}

    with open(fasta_library, "rt") as fasFile:
        current_seq = ""
        current_id = ""
        for line in fasFile:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    if current_seq:
                        seqs[current_id[:]] = current_seq[:]

                acc = acc_regex.search(line)
                if acc:
                    current_seq = ""
                    current_id = acc.group(0)
            else:
                current_seq += line
    input_file = r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\SCP-Data\LT-Hippo_1ug_SCP\Reports.pr_matrix.tsv"
    df = pd.read_csv(input_file, sep="\t")



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
