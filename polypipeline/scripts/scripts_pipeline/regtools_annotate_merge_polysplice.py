#!/usr/bin/env python
# coding: utf-8

######################################
# Deux paths à fournir:
## 'annotated_files_path' là où sont les fichiers annotés de regtools
## sur la dernière ligne, un path pour imprimer le tsv mergé
import pandas as pd
pd.set_option('display.max_columns', None)
import os

import getopt
import sys

opts, args = getopt.getopt(
    sys.argv[1:],
    'p:o:',
    ['path', 'outfile'],
)

for opt, arg in opts:
    if opt in ('-p', '--path'):
        annotated_files_path  = arg
    if opt in ('-o', '--outfile'):
        out_merge_file  = arg



#annotated_files_path = "/data-isilon/sequencing/ngs/NGS2024_7341/HG38_CNG/junctions/regtools/star/"

# Get list of files and their names from the path
# Replace the hyphens in the name by underscores

# Détecte le pattern de ma nomenclature de Regtools (un fichier annoté commence toujours par "annotated_")
pattern = "annotated_"
list_of_samples = []

# Récupère le path de chaque fichier et le nom de chaque patient
# Dans ce cas, je remplace aussi les '-' dans le nom du fichier par un '_'
# Je charges les fichiers et j'ajoute la colonne 'sample' qui contient le nom du patient (dans ce cas, je remets le nom sans les underscores sinon on a un message d'erreur quand on lance la mise en cache)
files_with_paths = [os.path.join(annotated_files_path, f) for f in os.listdir(annotated_files_path) if os.path.isfile(os.path.join(annotated_files_path, f))]
for file in files_with_paths:
    old_name = (os.path.splitext(os.path.basename(file))[0]).split('.')[0]
    if old_name.startswith(pattern):
        old_name = old_name[len(pattern):]
        file_name = old_name.replace("-", "_")
        globals()[file_name] = pd.read_csv(file, sep="\t")
        globals()[file_name]['sample'] = old_name
        list_of_samples.append(file_name)

# Je définie la fonction qui me sert à générer un dictionnaire de jonctions DA (on en a un par patient)
def DA_dict (df):
    DA = {}
    for index, row in df.iterrows():
        if row['anchor'] == 'DA':
            chr = row['chrom']
            start = row['start']
            end = row['end']
            strand = row['strand']
            score = row['score']
            name = row['name']
            if chr not in DA:
                DA[chr] = {'+': {'start': {}, 'end': {}}, '-': {'start': {}, 'end': {}}}
            if start not in DA[chr][strand]['start']:
                DA[chr][strand]['start'][start] = {}
            if end not in DA[chr][strand]['start'][start]:
                DA[chr][strand]['start'][start][end] = {'start': start, 'end': end, 'score': score, 'name': name}
            if end not in DA[chr][strand]['end']:
                DA[chr][strand]['end'][end] = {}
            if start not in DA[chr][strand]['end'][end]:
                DA[chr][strand]['end'][end][start] = {'start': start, 'end': end, 'score': score, 'name': name}
    return DA

# Je définie la fonction qui me permet de retrouver la jonction DA (si elle existe) pour chaque jonction D, A, NDA (pour le moment j'ignore les N dans cette analyse)
# Cette analyse dépend du anchor (D, A, NDA, N) puis du strand. En fonction de la combinaison anchor-strand, je sais si je dois chercher une correspondance en start ou en end et sur quel strand
def my_test(df, df_DA):
    df['reference_junction'] = None
    df['junc_normale_count'] = None
    df['DA_junction_s_'] = None
    df['ratio'] = None
    for index, row in df.iterrows():
        anchor = row['anchor']
        chr = row['chrom']
        start = row['start']
        end = row['end']
        strand = row['strand']
        score = row['score']
        if anchor == 'D' and strand == '+':
            if start in df_DA[chr][strand]['start']:
                max_end = 0
                max_score = 0
                list_DA_junc = []
                for key, value in df_DA[chr][strand]['start'][start].items():
                    list_DA_junc.append(df_DA[chr][strand]['start'][start][key]['name'])
                    if df_DA[chr][strand]['start'][start][key]['score'] > max_score:
                        max_score = df_DA[chr][strand]['start'][start][key]['score']
                        max_end = key
                df.at[index, 'DA_junction_s_'] = list_DA_junc
                df.at[index, 'reference_junction'] = df_DA[chr][strand]['start'][start][max_end]['name']
                df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['start'][start][max_end]['score']
                df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['start'][start][max_end]['score'])) * 100
            else:
                df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
                df.at[index, 'reference_junction'] = "No_matching_DA_junction"
                df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
                df.at[index, 'ratio'] = 100
        if anchor == 'D' and strand == '-':
            if end in df_DA[chr][strand]['end']:
                max_start = 0
                max_score = 0
                list_DA_junc = []
                for key, value in df_DA[chr][strand]['end'][end].items():
                    list_DA_junc.append(df_DA[chr][strand]['end'][end][key]['name'])
                    if df_DA[chr][strand]['end'][end][key]['score'] > max_score:
                        max_score = df_DA[chr][strand]['end'][end][key]['score']
                        max_start = key
                df.at[index, 'DA_junction_s_'] = list_DA_junc
                df.at[index, 'reference_junction'] = df_DA[chr][strand]['end'][end][max_start]['name']
                df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['end'][end][max_start]['score']
                df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['end'][end][max_start]['score'])) * 100
            else:
                df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
                df.at[index, 'reference_junction'] = "No_matching_DA_junction"
                df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
                df.at[index, 'ratio'] = 100
        if anchor == 'A' and strand == '+':
            if end in df_DA[chr][strand]['end']:
                max_start = 0
                max_score = 0
                list_DA_junc = []
                for key, value in df_DA[chr][strand]['end'][end].items():
                    list_DA_junc.append(df_DA[chr][strand]['end'][end][key]['name'])
                    if df_DA[chr][strand]['end'][end][key]['score'] > max_score:
                        max_score = df_DA[chr][strand]['end'][end][key]['score']
                        max_start = key
                df.at[index, 'DA_junction_s_'] = list_DA_junc
                df.at[index, 'reference_junction'] = df_DA[chr][strand]['end'][end][max_start]['name']
                df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['end'][end][max_start]['score']
                df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['end'][end][max_start]['score'])) * 100
            else:
                df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
                df.at[index, 'reference_junction'] = "No_matching_DA_junction"
                df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
                df.at[index, 'ratio'] = 100
        if anchor == 'A' and strand == '-':
            if start in df_DA[chr][strand]['start']:
                max_end = 0
                max_score = 0
                list_DA_junc = []
                for key, value in df_DA[chr][strand]['start'][start].items():
                    list_DA_junc.append(df_DA[chr][strand]['start'][start][key]['name'])
                    if df_DA[chr][strand]['start'][start][key]['score'] > max_score:
                        max_score = df_DA[chr][strand]['start'][start][key]['score']
                        max_end = key
                df.at[index, 'DA_junction_s_'] = list_DA_junc
                df.at[index, 'reference_junction'] = df_DA[chr][strand]['start'][start][max_end]['name']
                df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['start'][start][max_end]['score']
                df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['start'][start][max_end]['score'])) * 100
            else:
                df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
                df.at[index, 'reference_junction'] = "No_matching_DA_junction"
                df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
                df.at[index, 'ratio'] = 100
        if anchor == 'NDA':
            max_end = 0
            max_start = 0
            max_score_start = 0
            max_score_end = 0
            list_DA_junc = []
            if start in df_DA[chr][strand]['start']:
                for key, value in df_DA[chr][strand]['start'][start].items():
                    list_DA_junc.append(df_DA[chr][strand]['start'][start][key]['name'])
                    if df_DA[chr][strand]['start'][start][key]['score'] > max_score_start:
                        max_score_start = df_DA[chr][strand]['start'][start][key]['score']
                        max_end = key
            if end in df_DA[chr][strand]['end']:
                for key, value in df_DA[chr][strand]['end'][end].items():
                    list_DA_junc.append(df_DA[chr][strand]['end'][end][key]['name'])
                    if df_DA[chr][strand]['end'][end][key]['score'] > max_score_end:
                        max_score_end = df_DA[chr][strand]['end'][end][key]['score']
                        max_start = key
            if max_score_start > 0 or max_score_end > 0:
                if max_score_start >= max_score_end:
                    df.at[index, 'DA_junction_s_'] = list_DA_junc
                    df.at[index, 'reference_junction'] = df_DA[chr][strand]['start'][start][max_end]['name']
                    df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['start'][start][max_end]['score']
                    df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['start'][start][max_end]['score'])) * 100
                elif max_score_end > max_score_start:
                    df.at[index, 'DA_junction_s_'] = list_DA_junc
                    df.at[index, 'reference_junction'] = df_DA[chr][strand]['end'][end][max_start]['name']
                    df.at[index, 'junc_normale_count'] = df_DA[chr][strand]['end'][end][max_start]['score']
                    df.at[index, 'ratio'] = (score / (score + df_DA[chr][strand]['end'][end][max_start]['score'])) * 100
            else:
                df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
                df.at[index, 'reference_junction'] = "No_matching_DA_junction"
                df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
                df.at[index, 'ratio'] = 100
        if anchor == 'N':
            df.at[index, 'DA_junction_s_'] = "No_matching_DA_junction"
            df.at[index, 'reference_junction'] = "No_matching_DA_junction"
            df.at[index, 'junc_normale_count'] = "No_matching_DA_junction"
            df.at[index, 'ratio'] = 100
    return df

# Je loop dans ma liste de patients et je lance la fonction qui génère les dictionnaires de DA par patient. Je les renvoie dans le global context avec (globals[]) pour les réutiliser plus bas pour l'analyse.
list_of_samples_DA = []
for sample in list_of_samples:
    df = globals()[sample]
    df_DA_name = sample+"_DA"
    globals()[df_DA_name] = DA_dict(df)
    list_of_samples_DA.append(df_DA_name)

# Je loop dans ma liste de patients et je lance l'analyse. Les print c'est juste pour voir sur la console où j'en suis (comme c'est fait un par un), un équivalent de 'warn' dans ce cas.
for sample in list_of_samples:
    sample_DA = sample + "_DA"
    df = globals()[sample]
    df_DA = globals()[sample_DA]
    print(sample)
    print(sample_DA)
    my_test(df, df_DA)

# Je retire la colonne 'ratio' que j'utilisais sur ce projet pour fair du débugging, puis je réorganise les colonnes pour être dans le même ordre que les résults de Nico/ceux que tu parses déjà avec polysplice
for sample in list_of_samples:
    globals()[sample] = globals()[sample].drop(columns=['ratio'])
    globals()[sample] = globals()[sample][['name', 'chrom', 'start', 'end', 'score', 'strand', 'splice_site', 'acceptors_skipped', 'exons_skipped', 'donors_skipped', 'anchor', 'known_donor', 'known_acceptor', 'known_junction', 'gene_names', 'gene_ids', 'transcripts', 'DA_junction_s_', 'reference_junction', 'junc_normale_count', 'sample']]

# Je rédifinie le header pour ajouter le # en début
new_header = ["#name", "chr", "start", "end", "alt_count", "strand", "splice_site", "acceptors_skipped", "exons_skipped", "donors_skipped", "type", "known_donor", "known_acceptor", "known_junction", "gene", "ensid", "transcripts", "DA_junction_s_", "reference_junction", "junc_normale_count", "sample"]

# Je merge tous les patients pour générer un seul .tsv
list_of_dfs = [globals()[name] for name in list_of_samples]
test_df = pd.concat(list_of_dfs, ignore_index=True)
test_df.to_csv(out_merge_file, header= new_header, sep="\t", index = False)


