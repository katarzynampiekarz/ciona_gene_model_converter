import pandas as pd
import re
import numpy as np
from tkinter import *
from tkinter import messagebox
import pyperclip

"""
This program converts gene identifiers between KY and KH gene models for Ciona intestinalis.

The input can be a single gene or a batch input (one gene per line).
Example input format: KY21.Chr1.1. (just gene identifier closed with a dot, without transcription variant annotations).
KY21 genes go into the left input field, the right input field is meant for KH2013 genes.
The output is automatically copied to the clipboard.
Author: Katarzyna M Piekarz
Any issues can be reported to kpiekarz3@gatech.edu
Image source: FABA https://www.bpni.bio.keio.ac.jp/chordate/faba/1.4/top.html
Gene model source: http://ghost.zool.kyoto-u.ac.jp/default_ht.html
"""

#-----------------------------------------get_gene----------------------------------------------------------#

def clear():
    text_ky.delete('1.0', END)
    text_kh.delete('1.0', END)


def get_ky21_gene_batch(line, cutoff):
    kh2013 = pd.read_csv("kh2013_all.csv")
    ky21 = pd.read_csv("ky21_all.csv")

    gene_kh = line

    if gene_kh[0:2] != "KH":
        messagebox.showinfo(title="Not a proper KH2013 gene format", message="Provide a valid KH2013 gene name")
    else:
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
            text_ky.insert(END, eq)
        except:
            eq = "No match"
            text_ky.insert(END, eq)

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
                        eq = True
                        text_ky.insert(END, ky.group(0))
                        break

                if not eq:
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
                    text_ky.insert(END, eq)

            elif strands[0] == "-":
                ky21_chr = ky21_chr[ky21_chr['strand'] == "-"]
                kh2013_chr = kh2013_chr[kh2013_chr['strand'] == "-"]
                ends = []
                for index, row in kh2013_chr.iterrows():
                    if kh2013_chr['name'][index][0:n] == gene_kh:
                        ends.append(kh2013_chr['end'][index])

                for i in range(len(ends)):
                    if ends[i] in ky21_chr['end'].values:
                        rows = ky21_chr.loc[(ky21_chr['end'] == ends[i])]
                        ky_name = rows['name'].values[0]
                        ky = re.search(r".*?(?=v)", ky_name)
                        eq = True
                        text_ky.insert(END, ky.group(0))
                        break

                if not eq:
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

                        if abs(int(valR) - ends[i]) > abs(int(valL) - ends[i]):
                            val = valL
                        else:
                            val = valR

                        if abs(val - ends[i]) > cutoff:
                            eq = "No equivalent"

                        else:
                            rows = ky21_chr.loc[(ky21_chr['end'] == val)]
                            ky_name = rows['name'].values[0]
                            ky = re.search(r".*?(?=v)", ky_name)
                            ky_list_ends.append(ky.group(0))
                            eq = ky_list_ends[0]
                    text_ky.insert(END, eq)


def get_kh2013_gene_batch(line, cutoff):

    kh2013 = pd.read_csv("kh2013_all.csv")
    ky21 = pd.read_csv("ky21_all.csv")

    gene_ky21 = line

    if gene_ky21[0:4] != "KY21":
        messagebox.showinfo(title="Not a proper KY21 gene format", message="Provide a valid KY21 gene name")
    else:

        chrnum = re.search(r'KY21.(\w*).', gene_ky21).group(1)

        if kh2013['chrom'].str.contains(chrnum).any():
            ky21_chr = ky21[ky21['chrom'] == chrnum]
            kh2013_chr = kh2013[kh2013['chrom'] == chrnum]

            strands = []
            n = len(gene_ky21)
            eq = None

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
                        eq = True
                        text_kh.insert(END, kh.group(0))
                        break

                if not eq:
                    kh_sorted = np.sort(kh2013_chr['start'])
                    if len(kh_sorted) == 0:
                        eq = "No equivalent"
                    else:
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
                    text_kh.insert(END, eq)

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
                        eq = True
                        text_kh.insert(END, kh.group(0))
                        break

                if not eq:
                    kh_sorted = np.sort(kh2013_chr['end'])
                    if len(kh_sorted) == 0:
                        eq = "No equivalent"
                    else:
                        kh_list_ends = []
                        for i in range(len(ends)):
                            indR = np.searchsorted(kh_sorted, ends[i])
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
                    text_kh.insert(END, eq)

        else:
            text_kh.insert(END, "No equivalent")

def get_gene():
    cutoff = int(cutoff_entry.get())
    if text_kh.compare("end-1c", "==", "1.0"):
        genelist = text_ky.get("1.0", "end-1c").splitlines()
        for line in genelist:
            get_kh2013_gene_batch(line, cutoff)
            text_kh.insert(END, "\n")
            text_kh.update()
        pyperclip.copy(text_kh.get("1.0", END))

    elif text_ky.compare("end-1c", "==", "1.0"):
        genelist = text_kh.get("1.0", "end-1c").splitlines()
        for line in genelist:
            get_ky21_gene_batch(line, cutoff)
            text_ky.insert(END, "\n")
            text_kh.update()
        pyperclip.copy(text_ky.get("1.0", END))

#-----------------------------------------GUI---------------------------------------------------------------#

window = Tk()
window.title("Ciona Gene Model Converter | KY21 <-> KH2013")
window.config(padx=40, pady=40, bg="black")

canvas = Canvas(width=350, height=50, bg="black", highlightthickness=0)
logo_img = PhotoImage(file="ciona.png")
canvas.create_image(200, 20, image=logo_img)
canvas.grid(row=0, column=0, columnspan=3)

input_label_left = Label(text="KY21 example input format:\nKY21.Chr1.1.", font=("Arial", 9, "normal"), fg="white",
                         bg="black", justify=LEFT)
input_label_left.grid(column=0, row=1, columnspan=2)

input_label_right = Label(text="KH2013 example format:\nKH.C1.19.", font=("Arial", 9, "normal"), fg="white",
                          bg="black", justify=RIGHT)
input_label_right.grid(column=2, row=1)

generate_button = Button(text="Get equivalent gene", font=("Arial", 10, "normal"), width=18, command=get_gene)
generate_button.grid(column=0, row=3, columnspan=2, pady=10)

clear_button = Button(text="Clear input", font=("Arial", 10, "normal"), width=18, command=clear)
clear_button.grid(column=2, row=3)

text_ky = Text(height=5, width=19)
text_ky.focus()
text_ky.grid(column=0, row=2, pady=5, columnspan=2)

text_kh = Text(height=5, width=19)
text_kh.grid(column=2, row=2)

cutoff_label = Label(text="Distance cutoff:      ", font=("Arial", 9, "normal"), fg="white",
                         bg="black", justify=RIGHT)
cutoff_label.grid(column=0, row=4, sticky=E)

cutoff_entry = Entry(width=6, font=("Arial", 10, "normal"))
cutoff_entry.insert(END, string='500')
cutoff_entry.grid(column=1, row=4, sticky=W)


window.mainloop()
