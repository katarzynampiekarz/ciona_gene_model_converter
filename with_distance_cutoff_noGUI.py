import pandas as pd
import re
import numpy as np

"""
This program converts gene identifiers between KY and KH gene models for Ciona intestinalis.

The input file should be a csv file with a column containing genes to be converted, with column name either "kh" or "ky".
Example input format: KY21.Chr1.1. (just gene identifier closed with a dot, without transcription variant annotations).
Author: Katarzyna M Piekarz
Any issues can be reported to kpiekarz3@gatech.edu
Image source: FABA https://www.bpni.bio.keio.ac.jp/chordate/faba/1.4/top.html
Gene model source: http://ghost.zool.kyoto-u.ac.jp/default_ht.html
"""

# provide in parentheses below the path to your file with gene list
query = pd.read_csv("")
# you can adjust the distance cutoff (default: 500bp)
cutoff = 500


kh2013 = pd.read_csv("kh2013_all.csv")
ky21 = pd.read_csv("ky21_all.csv")


def get_ky21_gene_batch(gene_kh, cutoff):

    chrnums = []
    n = len(gene_kh)
    eq = None
    try:
        for index, row in kh2013.iterrows():
            if kh2013['name'][index][0:n] == gene_kh:
                chrnums.append(kh2013['chrom'][index])


        ky21_chr = ky21[ky21['chrom'] == chrnums[0]]
        kh2013_chr = kh2013[kh2013['chrom'] == chrnums[0]]
    except IndexError:
        eq = "No match in KH2013"
    except:
        eq = "No match"

    else:
        strands = []
        n = len(gene_kh)
        for index, row in kh2013_chr.iterrows():
            if kh2013_chr['name'][index][0:n] == gene_kh:
                strands.append(kh2013_chr['strand'][index])

        if strands[0] == "+":
            ky21_chr = ky21_chr[ky21_chr['strand'] == "+"]
            kh2013_chr = kh2013_chr[kh2013_chr['strand'] == "+"]
            starts = []
            for index, row in kh2013_chr.iterrows():
                if kh2013_chr['name'][index][0:n] == gene_kh:
                    starts.append(kh2013_chr['start'][index])

            for i in range(len(starts)):
                if starts[i] in ky21_chr['start'].values:
                    rows = ky21_chr.loc[(ky21_chr['start'] == starts[i])]
                    ky_name = rows['name'].values[0]
                    ky = re.search(r".*?(?=v)", ky_name)
                    eq = ky.group(0)
                    break

            if eq == None:
                ky_sorted = np.sort(ky21_chr['start'])
                ky_list = []
                for i in range(len(starts)):
                    indR = np.searchsorted(ky_sorted, starts[i])
                    if indR >= len(ky_sorted):
                        valR = ky_sorted[indR - 1]
                        valL = ky_sorted[indR - 2]
                    else:
                        valR = ky_sorted[indR]
                        valL = ky_sorted[indR-1]

                    if abs(int(valR) - starts[i]) > abs(int(valL) - starts[i]):
                        val = valL
                    else:
                        val = valR

                    if abs(int(val) - starts[i]) > cutoff:
                        eq = "No equivalent"

                    else:
                        rows = ky21_chr.loc[(ky21_chr['start'] == val)]
                        ky_name = rows['name'].values[0]
                        ky = re.search(r".*?(?=v)", ky_name)
                        ky_list.append(ky.group(0))

                eq = ky_list[0]
                #print(ky_list)

        elif strands[0] == "-":
            ky21_chr = ky21_chr[ky21_chr['strand'] == "-"]
            kh2013_chr = kh2013_chr[kh2013_chr['strand'] == "-"]
            ends = []
            for index, row in kh2013_chr.iterrows():
                if kh2013_chr['name'][index][0:n] == gene_kh:
                    ends.append(kh2013_chr['end'][index])
            #print("ends: ", ends)
            for i in range(len(ends)):
                if ends[i] in ky21_chr['end'].values:
                    rows = ky21_chr.loc[(ky21_chr['end'] == ends[i])]
                    ky_name = rows['name'].values[0]
                    ky = re.search(r".*?(?=v)", ky_name)
                    eq = ky.group(0)
                    break
            #print("eq: ", eq)
            if eq == None:
                ky_sorted = np.sort(ky21_chr['end'])
                ky_list_ends = []
                for i in range(len(ends)):
                    indR = np.searchsorted(ky_sorted, ends[i])
                    if indR >= len(ky_sorted):
                        valR = ky_sorted[indR - 1]
                        valL = ky_sorted[indR - 2]
                    else:
                        valR = ky_sorted[indR]
                        valL = ky_sorted[indR - 1]
                    #print("valR: ", valR, "valL: ", valL)
                    if abs(int(valR) - ends[i]) > abs(int(valL) - ends[i]):
                        val = valL
                    else:
                        val = valR
                    #print(val, ends[i])
                    if abs(val - ends[i]) > cutoff:
                        eq = "No equivalent"
                    else:
                        rows = ky21_chr.loc[(ky21_chr['end'] == val)]
                        ky_name = rows['name'].values[0]
                        ky = re.search(r".*?(?=v)", ky_name)
                        ky_list_ends.append(ky.group(0))
                    eq = ky_list_ends[0]
                    #print("ky_list_ends: ", ky_list_ends)
    finally:
        return eq

def get_kh2013_gene_batch(gene_ky21, cutoff):

    eq = None
    chrnum = re.search(r'KY21.(\w*).', gene_ky21).group(1)

    if kh2013['chrom'].str.contains(chrnum).any():
        ky21_chr = ky21[ky21['chrom'] == chrnum]
        kh2013_chr = kh2013[kh2013['chrom'] == chrnum]

        strands = []
        n = len(gene_ky21)

        for index, row in ky21_chr.iterrows():
            if ky21_chr['name'][index][0:n] == gene_ky21:
                strands.append(ky21_chr['strand'][index])

        if strands[0] == "+":
            kh2013_chr = kh2013_chr[kh2013_chr['strand'] == "+"]
            ky21_chr = ky21_chr[ky21_chr['strand'] == "+"]
            starts = []
            for index, row in ky21_chr.iterrows():
                if ky21_chr['name'][index][0:n] == gene_ky21:
                    starts.append(ky21_chr['start'][index])

            for i in range(len(starts)):
                if starts[i] in kh2013_chr['start'].values:
                    rows = kh2013_chr.loc[(kh2013_chr['start'] == starts[i])]
                    kh_name = rows['name'].values[0]
                    kh = re.search(r".*?(?=v)", kh_name)
                    eq = kh.group(0)
                    break

            if eq == None:
                kh_sorted = np.sort(kh2013_chr['start'])
                kh_list = []
                for i in range(len(starts)):
                    indR = np.searchsorted(kh_sorted, starts[i])
                    if indR >= len(kh_sorted):
                        valR = kh_sorted[indR - 1]
                        valL = kh_sorted[indR - 2]
                    else:
                        valR = kh_sorted[indR]
                        valL = kh_sorted[indR-1]

                    if abs(int(valR) - starts[i]) > abs(int(valL) - starts[i]):
                        val = valL
                    else:
                        val = valR

                    if abs(val - starts[i]) > cutoff:
                        eq = "No equivalent"

                    else:
                        rows = kh2013_chr.loc[(kh2013_chr['start'] == val)]
                        kh_name = rows['name'].values[0]
                        kh = re.search(r".*?(?=v)", kh_name)
                        kh_list.append(kh.group(0))
                        eq = kh_list[0]

        elif strands[0] == "-":
            kh2013_chr = kh2013_chr[kh2013_chr['strand'] == "-"]
            ky21_chr = ky21_chr[ky21_chr['strand'] == "-"]
            ends = []
            for index, row in ky21_chr.iterrows():
                if ky21_chr['name'][index][0:n] == gene_ky21:
                    ends.append(ky21_chr['end'][index])

            for i in range(len(ends)):
                if ends[i] in kh2013_chr['end'].values:
                    rows = kh2013_chr.loc[(kh2013_chr['end'] == ends[i])]
                    kh_name = rows['name'].values[0]
                    kh = re.search(r".*?(?=v)", kh_name)
                    eq = kh.group(0)
                    break

            if eq == None:
                kh_sorted = np.sort(kh2013_chr['end'])
                kh_list_ends = []
                for i in range(len(ends)):
                    indR = np.searchsorted(kh_sorted, ends[0])
                    if indR >= len(kh_sorted):
                        valR = kh_sorted[indR - 1]
                        valL = kh_sorted[indR - 2]
                    else:
                        valR = kh_sorted[indR]
                        valL = kh_sorted[indR - 1]

                    if abs(int(valR) - ends[i]) > abs(int(valL) - ends[i]):
                        val = valL
                    else:
                        val = valR

                    if abs(val - ends[i]) > cutoff:
                        eq = "No equivalent"

                    else:
                        rows = kh2013_chr.loc[(kh2013_chr['end'] == val)]
                        kh_name = rows['name'].values[0]
                        kh = re.search(r".*?(?=v)", kh_name)
                        kh_list_ends.append(kh.group(0))
                        eq = kh_list_ends[0]

    else:
        eq = "No equivalent"
    return eq

if query.columns[0] == "kh":
    genelist = query["kh"]
    input_model = "kh"
elif query.columns[0] == "ky":
    genelist = query["ky"]
    input_model = "ky"
else:
    print("Invalid input file")

genelist = genelist.tolist()

results = []
count = 1

for i in range(len(genelist)):
    gene_input = genelist[i]
    if input_model == "ky":
        gene_output = get_kh2013_gene_batch(gene_input, cutoff)
    elif input_model == "kh":
        gene_output = get_ky21_gene_batch(gene_input, cutoff)
    print(genelist[i], gene_output, count)
    results.append(gene_output)
    count += 1

if input_model == "kh":
    query["ky"] = results
else:
    query["kh"] = results

query.to_csv("converter_results.csv", index=False)
