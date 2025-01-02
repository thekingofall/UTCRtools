# utcrtools/options.py

import argparse
import configparser

def parse_args():
    parser = argparse.ArgumentParser(
        description="UTCRtools: TCR processing suite."
    )
    
    # Global arguments
    parser.add_argument("--FQ1", type=str, required=True, help="Fastq R1")
    parser.add_argument("--FQ2", type=str, required=True, help="Fastq R2")
    parser.add_argument("--BC", type=str, required=True, help="File with BC name and BC seq")
    parser.add_argument("--ini", type=str, required=True, help="config.ini")
    parser.add_argument("--FN", type=str, default="Filename", help="Filename prefix")
    parser.add_argument("--TCR", type=str, default="TRA,TRB", help="TCR types:<TRA|TRB|TRA,TRB|...>, default:TRA,TRB")
    parser.add_argument("--Num", type=str, default="3", help="TCR num, default:3")
    parser.add_argument("--BothSeq", type=str, default="CAGTGGTATCAACGCAGAG", help="Same sequence in your data, default:CAGTGGTATCAACGCAGAG")
    
    # Subparsers for different modules
    subparsers = parser.add_subparsers(dest='command', help='Available modules')
    subparsers.required = True  # Make subcommand required

    # Define each subcommand
    subparsers.add_parser('all', help="Transform + split BC + run Migec three times")
    subparsers.add_parser('transform', help="Only transform data")
    subparsers.add_parser('full', help="Transform + Migec + draw")
    subparsers.add_parser('full2', help="Transform + run Migec + draw")
    subparsers.add_parser('macaca', help="Run Migec for MacacaMulatta")
    subparsers.add_parser('macaca-transform', help="Transform + run Migec for MacacaMulatta")
    subparsers.add_parser('simple', help="Simple Migec run + draw")
    subparsers.add_parser('draw', help="Only draw")
    subparsers.add_parser('test', help="Test data")
    subparsers.add_parser('split-migec', help="Transform + BC_split + Migec")
    subparsers.add_parser('bc-migec-draw', help="BC_split + Migec + draw")
    subparsers.add_parser('migec-draw', help="Directly run Migec + draw on split files")
    subparsers.add_parser('run3-draw', help="Run Migec three times on split files + draw")
    subparsers.add_parser('migecfast', help="Migecfast + run Migec2 + draw")
    
    # Additional arguments for statistics
    parser.add_argument("--Statwhat", type=str, default="p", choices=['a', 'a1', 'p'], help="Data to stat: a (all), a1 (some), p (primary data)")
    parser.add_argument("--Statindex", type=str, help="Plot and stat folder. Run Statwhat module before if using this.")
    parser.add_argument("--DF", type=str, help="Folder to get statistics data")
    parser.add_argument("--SV", type=str, help="Folder to save statistics data")
    
    return parser.parse_args()

def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path, encoding='utf-8')
    software = dict(config.items("software"))
    print(software)
import configparser

def load_config(config_path):
    config = configparser.ConfigParser()
    config.optionxform = str  # 禁用键名的小写转换
    config.read(config_path, encoding='utf-8')
    software = dict(config.items("software"))
    print(software)
    return {
        "migecpath": software["migec"].strip('"'),
        "Drawpiepath": software["drawpie"].strip('"'),
        "PXTCR03_BCsplitpath": software["pxtcr03_bcsplit"].strip('"'),
        "PXTCR02_OSpath": software["pxtcr02_os"].strip('"'),
        "Draw_immu_index": software["draw_immu_index"].strip('"')
    }
