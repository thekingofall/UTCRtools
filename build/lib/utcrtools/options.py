# utcrtools/options.py

import argparse
import configparser

def parse_args():
    parser = argparse.ArgumentParser(description="UTCRtools: TCR processing suite.")
    
    # Filename arguments
    parser.add_argument("--FQ1", type=str, required=True, help="Fastq R1")
    parser.add_argument("--FQ2", type=str, required=True, help="Fastq R2")
    parser.add_argument("--BC", type=str, required=True, help="File with BC name and BC seq")
    parser.add_argument("--DF", type=str, help="Folder to get their stat")
    parser.add_argument("--SV", type=str, help="Folder to Save DF stat")
    parser.add_argument("--ini", type=str, required=True, help="config.ini")
    parser.add_argument("--FN", type=str, default="Filename", help="Filename you want to define")
    
    # TCR arguments
    parser.add_argument("--TCR", type=str, default="TRA,TRB", help="TCR types:<TRA|TRB|TRA,TRB|...>, default:TRA,TRB")
    parser.add_argument("--Num", type=str, default="3", help="TCR num, default:3")
    
    # Module arguments
    parser.add_argument("--Module", type=str, required=True, help="""Some Module you want to run:
a: transform + split BC + run migec in three times
f: transform only
fs: transform + migec + draw
fs2: transform + migecrun + draw
fsm: for MacacaMulatta
fsma: transform + MacacaMulatta
s: simple migec run + draw
t: only draw
o: testfordata
b: transform + BC_split + migec
bc: just BC_split + run migec + draw
bm: directly run migec + draw on splitted files
b3: run migecrun 3 times on splitted files + draw""")
    parser.add_argument("--Statwhat", type=str, default="p", help="""Some data you want to have a stat:
a: run all
a1: run some
p: only for primary data""")
    parser.add_argument("--BothSeq", type=str, default="CAGTGGTATCAACGCAGAG", help="""Some same Seq in your data:
default:CAGTGGTATCAACGCAGAG
Other seq you need""")
    parser.add_argument("--Statindex", type=str, help="""Plot: what folder you want to have a plot and stat,
but if want to use this module, you are better run the Statwhat module before""")
    
    return parser.parse_args()

def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path, encoding='utf-8')
    software = dict(config.items("software"))
    return {
        "migecpath": software["migec"].strip('"'),
        "Drawpiepath": software["Drawpie"].strip('"'),
        "PXTCR03_BCsplitpath": software["PXTCR03_BCsplit"].strip('"'),
        "PXTCR02_OSpath": software["PXTCR02_OS"].strip('"'),
        "Draw_immu_index": software["Draw_immu_index"].strip('"')
    }
