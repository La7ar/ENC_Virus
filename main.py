import pandas as pd
import numpy as np

def main():
    genome_enc = pd.read_csv("~/Downloads/genome_enc.csv")
    genes = pd.read_csv("~/Downloads/Gene_tsv.tsv", sep=",")
    genes.drop(genes.columns[0], axis=1, inplace=True)
    genes_with_complement = pd.read_csv("~/Downloads/Gene_tsv_with_complement.tsv", sep="\t")
    hmm = pd.read_csv("~/Downloads/RiboV1.4/RiboV1.4_HMMatches.tsv", sep="\t")
    # genes_genetic_codes = pd.read_csv("~/Downloads/gene_codes.csv")

    hmm_with_struct = classify_struct_of_genes(hmm)
    # genes['seq'] = genes_genetic_codes['x']
    complete_genes = remove_lines_of_uncompleted_genes(genes)

    # intersect(genes,hmms) - either struct or non-strcut (else --> "unknwon")
    # setdiff(genes,hmms) - "Unknown"
    # setdiff(hmm,genes) - -> /dev/null
    # setdiff(a,b) - in a a not in b


    hmm_with_struct['struct'].value_counts() #need for final report
    hmm_with_struct_without_unkown = hmm_with_struct[hmm_with_struct['struct'] != 'Unknown']


def count_number_of_unqiue_struct_columns_by_value(df):
    return df.groupby('struct')['struct'].count()
def classify_struct_of_genes(df):
    df2 = df.copy()
    list_of_struct = ["CP", "Coat", "Glycoprotein", "Nucleoprotein", "Membrane", "Capsid", "Envelope"]
    list_of_non_struct = ["RdRp", "Helicase", "Polymerase", "Methyltransferase", "Cap_MTase-GTase", "MTase", "GTase", "Pro", "Maturation"]
    df2['struct'] = np.where(df2['New_Name'].str.lower().str.contains('|'.join(map(str.lower, list_of_struct))), 'Structural',
                            np.where(df2['New_Name'].str.lower().str.contains('|'.join(map(str.lower, list_of_non_struct))), 'Non-Structural', 'Unknown'))
    return df2

def intersect_genes_by_ORFID(df1, df2):
    return df1[df1['ORFID'].isin(df2['ORFID'])]

def set_struct_column_to_unknown_if_ORFID_not_in_hmm(genes, hmm_with_struct):
    hmm_with_struct['struct'] = np.where(genes['ORFID'].isin(hmm_with_struct['ORFID']), hmm_with_struct['struct'], 'Unknown')
    return genes

def remove_row_if_ORFID_not_in_genes(hmm_with_struct, genes):
    return hmm_with_struct[hmm_with_struct['ORFID'].isin(genes['ORFID'])]

def remove_lines_of_uncompleted_genes(df):
    return df[~df['partial'].astype(str).isin(['01', '10', '1', '11'])]

if __name__ == '__main__':
    main()

