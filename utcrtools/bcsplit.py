import os
import sys
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def bcsplit(inputnamefile, inputfile):
    """
    根据 inputnamefile 中列出的 read ID，到 inputfile 中筛选对应的 read 并输出。
    """
    df_id = pd.read_csv(inputnamefile, names=["name"])
    name_list = set(df_id["name"].values)

    # 推断输出文件名
    FQ = inputfile.split(".")[1]  # e.g. "fastq"
    outputname = str(inputnamefile).split(".")[0]
    outname = outputname + "." + FQ + ".fastq"

    print(f"[bcsplit] Input name list: {inputnamefile}, input FASTQ: {inputfile}, output: {outname}")

    with open(inputfile) as f_in, open(outname, "w+") as out_handle:
        index = 0
        for title, seq, qual in FastqGeneralIterator(f_in):
            title1 = "@" + title.split(" ")[0]
            if title1 in name_list:
                seq_record = "\n".join(["@" + title, seq, "+", qual])
                out_handle.write(seq_record + "\n")
            index += 1

    print(f"[bcsplit] Processed {index} reads in file: {inputfile}")


def BCsplit_main():
    """

    usage: python bcsplit.py inputnamefile inputfile
    """
    if len(sys.argv) < 3:
        print("Usage: python bcsplit.py <inputnamefile> <inputfile>")
        sys.exit(1)

    inputnamefile = sys.argv[1]
    inputfile = sys.argv[2]
    bcsplit(inputnamefile, inputfile)
