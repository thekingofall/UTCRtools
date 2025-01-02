# utcrtools/choose.py

import os
import pandas as pd
import sys
import logging
import subprocess

from .utcrfunctions import (
    trans_data, BC_split, count_line, stat_module,
    transformer, migecrun, migecrun1, migecrun2, migecrun3,
    migecrunfast, migecfast, testfordata
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# 
def run_rscript(drawpiepath, input_file, output_file):
    run_command = ["Rscript", drawpiepath, input_file, output_file]
    logging.info(f"Running Rscript: {' '.join(run_command)}")
    try:
        subprocess.run(run_command, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Rscript failed with error: {e}")
        sys.exit(1)

# Define functions for each module
def run_all(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module all: transform + split BC + run Migec three times ...")
    trans_data(args=args, name=name)
    migecrun1(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    migecrun2(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)
    
    if args.BC:
        with open(args.BC, "r") as arr:
            for line in arr:
                fqname2, fqseq = line.strip().split("\t")
                for i in range(1, 4):
                    num2 = str(i)
                    migecrun(name=fqname2, num=num2, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
                    input_file2 = f"{fqname2}.cdrdata/{fqname2}_t{num2}TRA,TRB.cdrblast.txt"
                    output_file2 = f"{fqname2}.cdrdata/{fqname2}_{TCRtypes}_t{num2}"
                    run_rscript(cfg['Drawpiepath'], input_file2, output_file2)

def run_transform(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module transform: transform only ...")
    trans_data(args=args, name=name)

def run_full(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module full: transform + Migec + draw ...")
    trans_data(args=args, name=name)
    migecrun1(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    migecrun2(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_full2(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module full2: transform + run Migec + draw ...")
    trans_data(args=args, name=name)
    migecrun(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_macaca(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module macaca: run Migec for MacacaMulatta")
    migecrun3(name=name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta", migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/MacacaMulatta_{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/MacacaMulatta_{name}_TRB_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_macaca_transform(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module macaca-transform: transform + run Migec for MacacaMulatta")
    trans_data(args=args, name=name)
    migecrun1(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    migecrun3(name=name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta", migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/MacacaMulatta_{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/MacacaMulatta_{name}_TRB_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_simple(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module simple: simple Migec run + draw")
    migecrun1(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    migecrun2(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_draw(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module draw: only draw")
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_test(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module test: testfordata ...")
    testfordata(R=args.FQ1)

def run_split_migec(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module split-migec: transform + BC_split + Migec")
    trans_data(args=args, name=name)
    BC_split(BC=args.BC, fq1=fq1, fq2=fq2, PXTCR03_BCsplitpath=cfg["PXTCR03_BCsplitpath"])
    with open(args.BC, "r") as arr:
        for line in arr:
            fqname2, fqseq = line.strip().split("\t")
            migecrun1(name=fqname2, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
            migecrun2(name=fqname2, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
            input_file2 = f"{fqname2}.cdrdata/{fqname2}_t{num}TRA,TRB.cdrblast.txt"
            output_file2 = f"{fqname2}.cdrdata/{fqname2}_{TCRtypes}_t{num}"
            run_rscript(cfg['Drawpiepath'], input_file2, output_file2)

def run_bc_migec_draw(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module bc-migec-draw: BC_split + Migec + draw")
    BC_split(BC=args.BC, fq1=fq1, fq2=fq2, PXTCR03_BCsplitpath=cfg["PXTCR03_BCsplitpath"])
    with open(args.BC, "r") as arr:
        for line in arr:
            fqname2, fqseq = line.strip().split("\t")
            count_line(f"{fqname2}.FQ1.fastq")
            migecrun(name=fqname2, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
            input_file2 = f"{fqname2}.cdrdata/{fqname2}_t{num}TRA,TRB.cdrblast.txt"
            output_file2 = f"{fqname2}.cdrdata/{fqname2}_{TCRtypes}_t{num}"
            run_rscript(cfg['Drawpiepath'], input_file2, output_file2)

def run_migec_draw(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module migec-draw: directly run Migec + draw on split files")
    migecrun(name=name, num=num, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
    input_file = f"{name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt"
    output_file = f"{name}.cdrdata/{name}_{TCRtypes}_t{num}"
    run_rscript(cfg['Drawpiepath'], input_file, output_file)

def run_run3_draw(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module run3-draw: run Migec three times on split files + draw")
    with open(args.BC, "r") as arr:
        for line in arr:
            fqname2, fqseq = line.strip().split("\t")
            for i in range(1, 4):
                num2 = str(i)
                migecrun(name=fqname2, num=num2, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
                input_file2 = f"{fqname2}.cdrdata/{fqname2}_t{num2}TRA,TRB.cdrblast.txt"
                output_file2 = f"{fqname2}.cdrdata/{fqname2}_{TCRtypes}_t{num2}"
                run_rscript(cfg['Drawpiepath'], input_file2, output_file2)

def run_migecfast(args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext):
    logging.info("Module migecfast: migecfast + run Migec2 + draw")
    migecfast(BC2=args.BC, n=4, TCRtypes=TCRtypes, migecpath=cfg["migecpath"], PXTCR03_BCsplitpath=cfg["PXTCR03_BCsplitpath"])
    with open(args.BC, "r") as arr2:
        for line in arr2:
            fqname2, fqseq = line.strip().split("\t")
            for i in range(1, 4):
                num2 = str(i)
                migecrun2(name=fqname2, num=num2, TCRtypes=TCRtypes, migecpath=cfg["migecpath"])
                input_file2 = f"{fqname2}.cdrdata/{fqname2}_t{num2}TRA,TRB.cdrblast.txt"
                output_file2 = f"{fqname2}.cdrdata/{fqname2}_{TCRtypes}_t{num2}"
                run_rscript(cfg['Drawpiepath'], input_file2, output_file2)

def handle_statindex(args, cfg, name, num, fq1, global_Stattext):
    if args.Statindex == "a":
        logging.info("Module Statindex=a: stat the reads of fastq file")
        os.system(f"rm -rf {global_Stattext}")
        os.system(f"touch {global_Stattext}")
        data_tem = stat_module(FQone=args.FQ1, fq1=fq1, stat_file=global_Stattext, stat_what=args.Statwhat)
        logging.info(f'all is {data_tem[0]}, BC is {data_tem[1]}, Percentage is {data_tem[1]/data_tem[0]*100:.2f}%')

        if args.Statwhat == "a" and args.BC:
            with open(args.BC, "r") as arr, open(global_Stattext, "a+") as tem_data:
                for line in arr:
                    fqname2, fqseq = line.strip().split("\t")
                    fq_tem = count_line(f"{fqname2}.FQ1.fastq")[0]
                    tem_data.write(f"{fqname2},{fq_tem},{fq_tem/data_tem[1]*100:.2f}%\n")

    elif args.Statindex == "a1" and args.BC:
        logging.info("Module Statindex=a1: stat the Migec file")
        num3 = num
        dat2 = pd.DataFrame(columns=[
            "fqname_types","raw_types","raw_sum","clean_types","clean_sum",
            "clean_TRA_types","clean_TRA_sum","clean_TRB_types","clean_TRB_sum",
            "raw_max_cdr3_seq","clean_max_cdr3_seq",
            "raw_V_types","clean_V_types","clean_max_Vsegments","clean_max_V_num",
            "clean_J_types","clean_max_J_num","clean_max_Jsegments"
        ])

        with open(args.BC, "r") as arr:
            for line in arr:
                fqname2, fqseq = line.strip().split("\t")
                cdr_tem = pd.read_csv(f"{fqname2}.cdrdata/{fqname2}_t{num3}TRA,TRB.cdrblast.txt", sep="\t")
                logging.info(f'the types of {fqname2} raw is {cdr_tem.shape[0]}')
                cdr_clean = cdr_tem[~cdr_tem["CDR3 amino acid sequence"].str.contains("\?")]
                logging.info(f'the types of clean of {fqname2} is {cdr_clean.shape[0]}')
                cdr_clean.to_csv(f"{fqname2}.cdrdata/{fqname2}_t{num3}TRA,TRB.cdrblast.clean.csv", index=False)
                TRA = cdr_clean[cdr_clean["V segments"].str.contains("TRA")]
                logging.info(f'the types of TRA in {fqname2} is {TRA.shape[0]}')
                TRA.to_csv(f"{fqname2}.cdrdata/{fqname2}_t{num3}.TRA.clean.txt", sep="\t", index=False)

                TRB = cdr_clean[cdr_clean["V segments"].str.contains("TRB")]
                logging.info(f'the types of TRB in {fqname2} is {TRB.shape[0]}')
                TRB.to_csv(f"{fqname2}.cdrdata/{fqname2}_t{num3}.TRB.clean.txt", sep="\t", index=False)

                dat = {
                    "fqname_types": [fqname2],
                    "raw_types": [str(cdr_tem.shape[0])],
                    "raw_sum": [list(cdr_tem.iloc[:, [0]].sum())[0]],
                    "clean_types": [cdr_clean.shape[0]],
                    "clean_sum": [list(cdr_clean.iloc[:, [0]].sum())[0]],
                    "clean_TRA_types": [TRA.shape[0]],
                    "clean_TRA_sum": [list(TRA.iloc[:, [0]].sum())[0]],
                    "clean_TRB_types": [TRB.shape[0]],
                    "clean_TRB_sum": [list(TRB.iloc[:, [0]].sum())[0]],
                    "raw_max_cdr3_seq": cdr_tem.iloc[:, 3].iloc[0],
                    "clean_max_cdr3_seq": cdr_clean.iloc[:, 3].iloc[0],
                    "raw_V_types": [len(cdr_tem["V segments"].unique())],
                    "clean_V_types": [len(cdr_clean["V segments"].unique())],
                    "clean_max_Vsegments": [cdr_clean["V segments"].value_counts().idxmax()],
                    "clean_max_V_num": [cdr_clean["V segments"].value_counts().max()],
                    "clean_J_types": [len(cdr_clean["J segments"].unique())],
                    "clean_max_J_num": [cdr_clean["J segments"].value_counts().max()],
                    "clean_max_Jsegments": [cdr_clean["J segments"].value_counts().idxmax()]
                }
                dat = pd.DataFrame(dat)
                dat2 = dat2.append(dat.iloc[0, :], ignore_index=True)

        dat2.to_csv(f"{name}{num3}.BC_stat.csv", index=False)

    elif args.Statindex == "plot":
        folder_this = args.DF
        folder_this_name = args.SV
        run_rscript(cfg['Draw_immu_index'], folder_this, folder_this_name)
