import io
import matplotlib
import sys
import sklearn.neighbors
import time
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, normalize, RobustScaler

sys.modules['sklearn.neighbors.base'] = sklearn.neighbors._base
matplotlib.use('agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.ioff()

import pandas as pd
import pandera as pa
import os
import re
import numpy as np
from glob import glob

from scipy.stats import variation, sem
import scipy.cluster.hierarchy as sch

from unimod_mapper import UnimodMapper

from sequal.sequence import Sequence

from uniprotparser.parser import acc_regex, UniprotParser
from pandas import Series
from plotnine import ggplot, aes, geom_col, scale_fill_brewer, facet_wrap, theme_minimal, geom_text, geom_boxplot, \
    theme, element_text, save_as_pdf_pages, ggtitle, geom_density, scale_color_brewer, geom_line, geom_point, \
    geom_errorbar, geom_tile

import seaborn as sns
from missingpy import MissForest

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly
    correlated variables are next to eachother

    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max() / 2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold,
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)

    if not inplace:
        corr_array = corr_array.copy()

    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]


def skip_rows(x):
    row = x.strip()
    if row == "":
        return False
    if row.starswith("---"):
        return False



def initial_check_pr(df: pd.DataFrame, samples: list[str]):
    schema = {
        "Protein.Group": pa.Column(str),
        "Protein.Ids": pa.Column(str),
        "Protein.Names": pa.Column(str, nullable=True),
        "Genes": pa.Column(str, nullable=True),
        #"First.Protein.Description": pa.Column(str, nullable=True),
        "Proteotypic": pa.Column(int),
        "Stripped.Sequence": pa.Column(str),
        "Modified.Sequence": pa.Column(str),
        "Precursor.Charge": pa.Column(int),
        "Precursor.Id": pa.Column(str)
    }
    for c in samples:
        schema[c] = pa.Column(float, nullable=True)

    pa_schem = pa.DataFrameSchema(schema)
    pa_schem(df)


def initial_check_pg(df: pd.DataFrame, samples: list[str]):
    schema = {
        "Protein.Group": pa.Column(str),
        "Protein.Ids": pa.Column(str),
        "Protein.Names": pa.Column(str, nullable=True),
        "Genes": pa.Column(str, nullable=True),
        #"First.Protein.Description": pa.Column(str, nullable=True),
    }
    for c in samples:
        schema[c] = pa.Column(float, nullable=True)

    pa_schem = pa.DataFrameSchema(schema)
    pa_schem(df)


def get_count_null(df: pd.DataFrame, samples: list[str]):
    data = []
    for s in samples:
        s_na = df[df[s].notnull()]
        data.append([s, s_na.shape[0], "Not Null"])
        data.append([s, df.shape[0] - s_na.shape[0], "Null"])
    na_df = pd.DataFrame(data, columns=["Sample", "Count", "Category"])
    na_df["Category"] = pd.Categorical(na_df["Category"].values, ordered=True, categories=["Not Null", "Null"])
    return na_df


class File:
    def __init__(self, name):
        self.name = name
        self.alias = name
        self.data = pd.DataFrame()

    def set_data(self, data):
        self.data = data


class Diann:
    def __init__(self, folder_path, temp_folder_path=".", fasta_lib_path="", goa_file="", go_obo=""):
        self.temp_folder_path = temp_folder_path
        os.makedirs(self.temp_folder_path, exist_ok=True)
        self.progress_file = open(os.path.join(self.temp_folder_path, "progress.txt"), "wt")
        self.folder_path = folder_path
        self.fasta_lib_path = fasta_lib_path
        self.goa_file = goa_file
        self.go_obo = go_obo
        self.annotation = {}
        self.raw_file_location = self.parse_log()
        self.pr = {}
        self.pg = {}
        self.coverage = pd.DataFrame()
        self.modification = pd.DataFrame()

        self.process_data()
        self.joined_pr = self.join_df(self.pr)
        self.joined_pg = self.join_df(self.pg)
        self.unimod_mapper = UnimodMapper()
        if self.fasta_lib_path == "":
            self.fasta_lib = self.get_fasta_lib()

        else:
            self.fasta_lib = self.load_fasta_library()
        self.draw_null()
        self.pca(self.joined_pg, self.temp_folder_path)
        self.draw_profile(self.joined_pg, os.path.join(self.temp_folder_path, "profile.pg.pdf"))
        self.draw_total_intensity(self.joined_pr, os.path.join(self.temp_folder_path, "total.intensity.pr.pdf"))
        self.draw_detected_genes(self.joined_pr, os.path.join(self.temp_folder_path, "gene.detected.pr.pdf"))
        self.draw_cv(self.joined_pr, os.path.join(self.temp_folder_path, "cv.pr.pdf"),
                     ["Condition", "Protein.Group", "Modified.Sequence"])
        self.draw_cv(self.joined_pg, os.path.join(self.temp_folder_path, "cv.pg.pdf"), ["Condition", "Protein.Group"])
        self.protein_coverage()
        self.modification_map()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.progress_file.close()

    def write_progress(self, message):
        self.progress_file.write(str(time.time()) + "\t" + message + "\n")
        self.progress_file.flush()

    def join_df(self, df_dict):
        combined_pr = []
        order = []
        condition_order = []
        for k in df_dict:
            if k in self.annotation:
                folder = k
                df = pd.melt(df_dict[k],
                             id_vars=[i for i in df_dict[k].columns if i not in self.raw_file_location[k]],
                             value_vars=self.raw_file_location[k],
                             value_name="Intensity",
                             var_name="Sample")
                label = []
                for i, r in df.iterrows():
                    if r["Sample"] in self.annotation[folder]:
                        df.at[i, "Condition"] = self.annotation[folder][r["Sample"]]
                        label.append(self.annotation[folder][r["Sample"]])
                        if self.annotation[folder][r["Sample"]] not in condition_order:
                            condition_order.append(self.annotation[folder][r["Sample"]])
                combined_pr.append(df)
            else:
                folder, condition = os.path.split(k)
                condition_order.append(condition)
                df_dict[k]["Condition"] = Series([condition] * len(df_dict[k].index), index=df_dict[k].index)

                df = pd.melt(df_dict[k],
                             id_vars=[i for i in df_dict[k].columns if i not in self.raw_file_location[k]],
                             value_vars=self.raw_file_location[k],
                             value_name="Intensity",
                             var_name="Sample")
                df["Location"] = k
                combined_pr.append(df)
            order = order + self.raw_file_location[k]
        if len(combined_pr) < 2:
            df = combined_pr[0]
        else:
            df = pd.concat(combined_pr, ignore_index=True)
        df["Sample"] = pd.Categorical(df["Sample"].values, categories=order)
        df["Condition"] = pd.Categorical(df["Condition"].values, categories=condition_order)
        return df

    def parse_log(self):
        file_location = {}
        self.write_progress("Began parsing log file")
        for l in glob(os.path.join(self.folder_path, "**", "*"), recursive=True):
            if l.endswith("Reports.log.txt"):
                folder, _ = os.path.split(l)
                file_location[folder] = []
                annotation_path = os.path.join(folder, "annotation.txt")
                if os.path.exists(annotation_path):
                    self.annotation[folder] = {}
                    df = pd.read_csv(annotation_path, sep="\t")
                    for i, r in df.iterrows():
                        self.annotation[folder][r["Sample"]] = r["Condition"]
                    file_location[folder] = list(df["Sample"].values)

                else:
                    with open(l, "rt") as log_stuff:
                        for line in log_stuff:
                            line = line.strip()
                            if "Loading run" in line:
                                match = re.search("Loading run (.+)", line)
                                if match:
                                    file = match.group(1)
                                    if file not in file_location[folder]:
                                        file_location[folder].append(file)
        self.write_progress("Completed parsing log file")
        return file_location

    def process_pr(self, folder_path, file_list):
        self.write_progress("Began processing pr_matrix file")
        pr_file = os.path.join(folder_path, "Reports.pr_matrix.tsv")
        df = pd.read_csv(pr_file, sep="\t")
        initial_check_pr(df, file_list)
        self.write_progress("Completed processing pr_matrix file")
        return df

    def process_pg(self, folder_path, file_list):
        self.write_progress("Began processing pg_matrix file")
        pg_file = os.path.join(folder_path, "Reports.pg_matrix.tsv")
        df = pd.read_csv(pg_file, sep="\t")
        initial_check_pg(df, file_list)
        print(self.goa_file, self.go_obo)
        if self.goa_file != "" and self.go_obo != "":
            self.write_progress("Began adding gene ontology definitions to pg_matrix file")
            a = set()
            dic = {}
            for prs in df["Protein.Group"]:
                dic[prs] = prs.split(";")
                for pr in prs.split(";"):
                    a.add(pr)

            g = GafReader(self.goa_file)
            gobo = GODag(self.go_obo)

            annotation = g.get_ns2assc()
            gos = {"BP": [], "CC": [], "MF": []}

            for prs in df["Protein.Group"]:
                go = {"BP": set(), "CC": set(), "MF": set()}
                for p in dic[prs]:
                    for ns in annotation:
                        if p in annotation[ns]:
                            for i in annotation[ns][p]:
                                if i in gobo:
                                    go[ns].add(gobo[i].name)
                for ns in go:
                    gos[ns].append(";".join(go[ns]))

            for ns in gos:
                df["GO.{}".format(ns)] = Series(gos[ns], index=df.index)
            self.write_progress("Completed adding gene ontology definitions to pg_matrix file")
        self.write_progress("Completed processing pg_matrix file")
        return df

    def process_data(self):
        for l in self.raw_file_location:
            self.pr[l] = self.process_pr(l, self.raw_file_location[l])
            self.pg[l] = self.process_pg(l, self.raw_file_location[l])
            _, folder = os.path.split(l)
            os.makedirs(os.path.join(self.temp_folder_path, folder), exist_ok=True)

            self.draw_correlation_matrix(self.pr[l], os.path.join(self.temp_folder_path, folder))
            peptide = self.pr[l].groupby(["Protein.Group"])["Proteotypic"].size()
            peptide = peptide.reset_index()
            peptide = peptide.rename(columns={peptide.columns[-1]: "Peptide.Count"})
            unique_peptide = self.pr[l][self.pr[l]["Proteotypic"] == 1].groupby(["Protein.Group"]).size()
            unique_peptide = unique_peptide.reset_index()
            unique_peptide = unique_peptide.rename(columns={unique_peptide.columns[-1]: "Unique.Peptide.Count"})
            unique_peptide = unique_peptide
            peptide = peptide.merge(unique_peptide, on="Protein.Group")
            self.pg[l].merge(peptide, on="Protein.Group").to_csv(
                os.path.join(self.temp_folder_path, folder, "Reports.pg_matrix.with.peptide.count.tsv"), sep="\t", index=False)

    def draw_null(self):
        for l in self.raw_file_location:
            self.write_progress("Began checking for null values")
            _, folder = os.path.split(l)
            os.makedirs(os.path.join(self.temp_folder_path, folder), exist_ok=True)
            self.draw_report_stats(l, os.path.join(self.temp_folder_path, folder))
            folder = os.path.join(self.temp_folder_path, folder)
            self._draw_null(
                l,
                data=self.pr[l],
                path=os.path.join(folder, "temp.pr_matrix.countnull.pdf"),
                file_list=self.raw_file_location[l])
            self._draw_null(
                l,
                data=self.pg[l],
                path=os.path.join(folder, "temp.pg_matrix.countnull.pdf"),
                file_list=self.raw_file_location[l])
            self.write_progress("Completed checking for null values")

    def draw_report_stats(self, folder, path):
        self.write_progress("Began drawing plot for protein identified number")
        report_path = os.path.join(folder, "Reports.stats.tsv")
        df = pd.read_csv(report_path, sep="\t")
        for i, r in df.iterrows():
            sample = self.raw_file_location[folder][i]
            df.at[i, "Condition"] = self.annotation[folder][sample]
        mean_data = []
        for label, group_data in df.groupby("Condition"):
            mean_data.append(
                [label, np.mean(group_data['Proteins.Identified']), sem(group_data["Proteins.Identified"])])
        df_mean = pd.DataFrame(mean_data, columns=["Condition", "Count", "Error"])
        result = ggplot() + \
                 geom_col(df_mean, aes(x="Condition", y="Count", fill="Condition")) + \
                 geom_point(df, aes(x="Condition", y="Proteins.Identified"), fill="black", position="jitter") + \
                 geom_errorbar(
                     df_mean, aes(x="Condition", ymin="Count-Error", ymax="Count+Error")) + \
                 scale_fill_brewer(type="qual", palette="Pastel1") + theme_minimal()
        result.save(os.path.join(path, "proteins.identified.pdf"))
        self.write_progress("Completed drawing plot for protein identified number")

    def _draw_null(self, folder, data, path, file_list):
        data_pr = get_count_null(data, file_list)
        plots = []
        mean_data = []
        for label, group_data in data_pr.groupby("Sample"):
            plots.append(
                ggplot(group_data, aes(x='Category', y='Count', fill='Category', label='Count')) \
                + scale_fill_brewer(type="qual", palette="Pastel2") \
                + geom_col() \
                + geom_text() \
                + theme_minimal() + theme(figure_size=(8, 6)) + ggtitle(label))
        not_null = data_pr[data_pr["Category"] == "Not Null"]
        not_null["Condition"] = Series([self.annotation[folder][i] for i in not_null["Sample"]], index=not_null.index)
        for label, group_data in not_null.groupby("Condition"):
            mean_data.append([label, np.mean(group_data['Count']), sem(group_data["Count"])])
        df_mean = pd.DataFrame(mean_data, columns=["Condition", "Count", "Error"])
        plots = [
                    ggplot() +
                    geom_col(df_mean, aes(x="Condition", y="Count", fill="Condition")) +
                    geom_point(not_null, aes(x="Condition", y="Count"), fill="black", position="jitter") +
                    geom_errorbar(
                        df_mean, aes(x="Condition", ymin="Count-Error", ymax="Count+Error")) +
                    scale_fill_brewer(type="qual", palette="Pastel1") + theme_minimal()] + plots
        save_as_pdf_pages(plots, path)
        # result.save(path, "pdf", dpi=600, width=7, height=10)

    def draw_profile(self, data, path):
        self.write_progress("Began drawing profile plot")
        d = data.copy()
        d["log2Intensity"] = np.log2(d["Intensity"])
        result = ggplot(d, aes(x='Sample', y='log2Intensity', fill='Condition')) \
                 + geom_boxplot() + theme_minimal() + theme(axis_text_x=element_text(rotation=90))
        result.save(path, "pdf", width=10, height=15)
        self.write_progress("Completed drawing profile plot")

    def draw_total_intensity(self, data, path):
        self.write_progress("Began drawing total intensity plot")
        total_intensity = data.groupby(["Sample", "Condition"])["Intensity"].sum()
        total_intensity = total_intensity.reset_index()
        result = ggplot(total_intensity, aes(x="Sample", y="Intensity", fill="Condition")) \
                 + geom_col(colour="black") + theme_minimal() + theme(axis_text_x=element_text(rotation=90))
        result.save(path, "pdf", width=10, height=15)
        self.write_progress("Completed drawing total intensity plot")

    def draw_detected_genes(self, data, path):
        genes = data[data["Intensity"].notnull()]
        genes = genes.groupby(["Sample", "Condition"])["Genes"].nunique()
        genes = genes.reset_index()
        result = ggplot(genes, aes(x="Sample", y="Genes", fill="Condition")) \
                 + geom_col(colour="black") + theme_minimal() + theme(axis_text_x=element_text(rotation=90))
        result.save(path, "pdf", width=10, height=15)

    def load_fasta_library(self):
        seqs = {}

        with open(self.fasta_lib_path, "rt") as fasFile:
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

            seqs[current_id[:]] = current_seq[:]
        return seqs

    def get_seq(self, proteins):
        for p in proteins["Protein.Group"].split(";"):
            if p in self.fasta_lib:
                if p in self.fasta_lib:
                    yield p, self.fasta_lib[p], proteins["Genes"]
                else:
                    yield p, None, proteins["Genes"]

    def protein_coverage(self):
        self.write_progress("Began analyzing protein coverage")
        cover_dict = {}
        gene_dict = {}
        for i, r in self.joined_pr.iterrows():
            if pd.notnull(r["Intensity"]):
                for protein, seq, gene in self.get_seq(r):
                    gene_dict[protein] = gene
                    if protein not in cover_dict:
                        cover_dict[protein] = {}
                    if r["Sample"] not in cover_dict[protein]:
                        cover_dict[protein][r["Sample"]] = {}
                    if seq:
                        try:
                            ind = seq.index(r["Stripped.Sequence"])
                        except ValueError:
                            if protein+"-1" in self.fasta_lib:
                                seq = self.fasta_lib[protein+"-1"]
                                ind = seq.index(r["Stripped.Sequence"])
                            else:
                                ind = -1
                        if ind > -1:
                            for n in range(len(r["Stripped.Sequence"])):
                                pos = ind + n
                                if pos not in cover_dict[protein][r["Sample"]]:
                                    cover_dict[protein][r["Sample"]][pos] = 0
                                cover_dict[protein][r["Sample"]][pos] = cover_dict[protein][r["Sample"]][pos] + 1
        result = []
        for p in cover_dict:
            if p in self.fasta_lib:
                for s in cover_dict[p]:
                    covered_pos = list(cover_dict[p][s].keys())
                    seq_len = len(self.fasta_lib[p])
                    result.append([p, len(covered_pos), seq_len, len(covered_pos) / seq_len, gene_dict[p], s])

        df = pd.DataFrame(result,
                          columns=["Protein", "Length Covered", "Length", "Covered Proportion", "Gene names", "Sample"])
        self.coverage = df
        self.coverage.to_csv(os.path.join(self.temp_folder_path, "coverage.txt"), sep="\t", index=False)
        self.write_progress("Completed analyzing protein coverage")

    def modification_map(self):
        self.write_progress("Began mapping modification sites")
        mod_dict = {}
        result = []
        for i, r in self.joined_pr.iterrows():
            if pd.notnull(r["Intensity"]):
                for protein, seq, gene in self.get_seq(r):

                    if protein not in mod_dict:
                        mod_dict[protein] = {}
                    if r["Sample"] not in mod_dict[protein]:
                        mod_dict[protein][r["Sample"]] = {}

                    if seq:
                        try:
                            ind = seq.index(r["Stripped.Sequence"])
                        except ValueError:
                            if protein+"-1" in self.fasta_lib:
                                seq = self.fasta_lib[protein+"-1"]
                                ind = seq.index(r["Stripped.Sequence"])
                            else:
                                ind = -1
                        modded_pep = Sequence(r["Modified.Sequence"])
                        for n, a in enumerate(modded_pep):
                            if len(a.mods) > 0:
                                for m in a.mods:
                                    name = self.unimod_mapper.id_to_name(str(m).replace("UniMod:", ""))
                                    pos = ind + n
                                    result.append([protein, pos, name[0], gene, r["Sample"], r["Intensity"]])

        df = pd.DataFrame(result, columns=["Protein", "Position", "Modification", "Gene names", "Sample", "Intensity"])
        res = df.groupby(["Protein", "Position", "Modification", "Gene names", "Sample"]).size()
        total_intensity = df.groupby(["Protein", "Position", "Modification", "Gene names", "Sample"])["Intensity"].sum()
        res = pd.concat([res, total_intensity], axis=1)
        self.modification = res.reset_index()
        intensity_matrix = self.modification[
            ["Protein", "Position", "Modification", "Gene names", "Sample", "Intensity"]]
        intensity_matrix.set_index(["Protein", "Position", "Modification", "Gene names", "Sample"], inplace=True)
        intensity_matrix = intensity_matrix.unstack("Sample")
        intensity_matrix.to_csv(os.path.join(self.temp_folder_path, "modification.map.intensity.matrix.txt"), sep="\t",
                                index=True)
        self.modification.rename(columns={0: "Instance counts"}, inplace=True)
        self.modification.to_csv(os.path.join(self.temp_folder_path, "modification.map.txt"), sep="\t", index=False)
        self.write_progress("Completed mapping modification sites")

    def calculate_CV(self, joined_df: pd.DataFrame, group):
        arr = []
        for i, g in joined_df.groupby(group):
            a = variation(g["Intensity"], nan_policy="omit")
            arr.append([i[0], a])
        df = pd.DataFrame(arr, columns=["Condition", "CV"])
        df = df[pd.notnull(df["CV"]) & (df["CV"] != "--")]
        return df

    def draw_cv_density(self, df: pd.DataFrame, filename):
        self.write_progress("Began creating CV plot")
        plots = []
        custom_legend = {}
        with PdfPages(filename) as pdf_pages:

            for label, group_data in df.groupby(["Condition"]):
                fig, ax = plt.subplots()
                title = label + " ({})".format(np.round(group_data["CV"].median(), 2))
                custom_legend[label] = title
                sns.kdeplot(data=group_data, x="CV", linewidth=5, ax=ax).set(title=title)
                plots.append(fig)
            fig, ax = plt.subplots()
            sns.kdeplot(data=df.assign(Label=df["Condition"].map(custom_legend)), hue="Label", x="CV", linewidth=5,
                        ax=ax)
            plots = [fig] + plots
            for p in plots:
                pdf_pages.savefig(p)
        # for label, group_data in df.groupby(["Condition"]):
        #     title = label + " ({})".format(np.round(group_data["CV"].median(), 2))
        #     custom_legend.append(title)
        #     plots.append(
        #         ggplot(group_data, aes(x="CV"))
        #         + geom_line(size=2, stat="density")
        #         # + geom_density(alpha=0.1)
        #         + theme_minimal() + ggtitle(title))
        #
        # plots = [ggplot(df, aes(x="CV", colour="Condition"))
        #          + geom_line(size=2, stat="density")
        #          # + geom_density(alpha=0.1)
        #          + scale_color_brewer(
        #     type="qual", palette="Pastel1", labels=custom_legend) + theme_minimal()] + plots
        # save_as_pdf_pages(plots, filename)
        self.write_progress("Completed creating CV plot")

    def draw_cv(self, joined_df: pd.DataFrame, filename, group):
        cv = self.calculate_CV(joined_df, group)
        cv.to_csv(filename + ".txt", sep="\t", index=False)
        self.draw_cv_density(cv, filename)

    def draw_correlation_matrix(self, df: pd.DataFrame, path):
        self.write_progress("Began creating correlation matrix")
        l, folder = os.path.split(path)
        corr_mat = df[df.columns[10:]].corr()
        without_cluster_mat = corr_mat.reset_index()
        without_cluster_mat = without_cluster_mat.rename(columns={without_cluster_mat.columns[0]: "X"})
        orderX = [i for i in without_cluster_mat["X"]]
        orderY = [i for i in without_cluster_mat.columns if i != "X"]
        melted = pd.melt(without_cluster_mat,
                         id_vars=["X"],
                         value_vars=[i for i in without_cluster_mat.columns if i != "X"],
                         value_name="Intensity",
                         var_name="Y"
                         )
        melted["X"] = pd.Categorical(melted["X"], categories=orderX)
        melted["Y"] = pd.Categorical(melted["Y"], categories=orderY)
        plots = ggplot(melted, aes(x="X", y="Y", fill="Intensity")) + geom_tile() + theme_minimal() \
                + theme(axis_text_x=element_text(rotation=90, size=5))
        plots.save(os.path.join(path, "correlation_matrix.without_cluster.pdf"))

        clustered = cluster_corr(corr_mat)
        clustered.to_csv(os.path.join(path, "correlation_matrix.tsv"), sep="\t")
        clustered = clustered.reset_index()
        clustered = clustered.rename(columns={clustered.columns[0]: "X"})
        orderX = [i for i in clustered["X"]]
        orderY = [i for i in clustered.columns if i != "X"]
        melted = pd.melt(clustered,
                         id_vars=["X"],
                         value_vars=[i for i in clustered.columns if i != "X"],
                         value_name="Intensity",
                         var_name="Y"
                         )
        melted["X"] = pd.Categorical(melted["X"], categories=orderX)
        melted["Y"] = pd.Categorical(melted["Y"], categories=orderY)
        plots = ggplot(melted, aes(x="X", y="Y", fill="Intensity")) + geom_tile() + theme_minimal() \
                + theme(axis_text_x=element_text(rotation=90, size=5))
        plots.save(os.path.join(path, "correlation_matrix.pdf"))
        self.write_progress("Completed creating correlation matrix")

    def get_fasta_lib(self):
        self.write_progress("Began getting sequence data from UniProt database")
        accessions = set()
        for i in self.joined_pg["Protein.Group"]:
            accs = i.split(";")
            for a in accs:
                accessions.add(a)

        p = UniprotParser(accessions)
        seqs = {}
        for res in p.parse():
            current_seq = ""
            current_id = ""
            for line in io.StringIO(res):
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
            seqs[current_id[:]] = current_seq[:]
        self.write_progress("Completed getting sequence data from UniProt database")
        return seqs

    def pca(self, df: pd.DataFrame, path):
        self.write_progress("Began PCA analysis")
        d = df.copy()
        d["Sample"] = d["Sample"].astype(str)
        a = d["Sample"].unique()
        d = d[[c for c in d.columns if c != "Condition"]]
        d = d.pivot(columns="Sample", values="Intensity", index=[c for c in d.columns if c not in ["Sample", "Intensity"]])
        imputer = MissForest()
        imputed = imputer.fit_transform(d[a])
        scaler = RobustScaler()
        imputed = scaler.fit_transform(imputed)
        pca = PCA(n_components=2)
        pc = pca.fit_transform(imputed.T)
        pc1 = "PC1({})".format(round(pca.explained_variance_ratio_[0], 2))
        pc2 = "PC2({})".format(round(pca.explained_variance_ratio_[1], 2))
        new_df = pd.DataFrame(data=pc, columns=[pc1, pc2])
        new_df["Sample"] = pd.Series(a, index=new_df.index)
        for i, r in new_df.iterrows():
            for l in self.annotation:
                if r["Sample"] in self.annotation[l]:
                    new_df.at[i, "Condition"] = self.annotation[l][r["Sample"]]
                    break
        new_df.to_csv(os.path.join(path, "pca.csv"), index=False)
        plot = ggplot(new_df, aes(x=pc1, y=pc2, fill="Condition", colour="Condition")) + geom_point(alpha=0.50, size=2) + theme_minimal()
        plot.save(os.path.join(path, "pca.pdf"))
        self.write_progress("Completed PCA analysis")
