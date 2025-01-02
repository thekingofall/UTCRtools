# utcrtools/main.py

import os
import sys
import datetime
import logging
from .options import parse_args, load_config
from .utcrchoose import (
    run_all, run_transform, run_full, run_full2, run_macaca,
    run_macaca_transform, run_simple, run_draw, run_test,
    run_split_migec, run_bc_migec_draw, run_migec_draw,
    run_run3_draw, run_migecfast, handle_statindex
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # Parse arguments
    args = parse_args()

    starttime = datetime.datetime.now()

    # Load configuration
    cfg = load_config(args.ini)
    migecpath = cfg["migecpath"]
    Drawpiepath = cfg["Drawpiepath"]
    PXTCR03_BCsplitpath = cfg["PXTCR03_BCsplitpath"]
    PXTCR02_OSpath = cfg["PXTCR02_OSpath"]
    Draw_immu_index = cfg["Draw_immu_index"]

    # Initialize variables
    R = args.FQ1
    R2 = args.FQ2
    TCRtypes = args.TCR
    num = str(args.Num)
    name = args.FN  # If not provided, default is "Filename"
    Sameseq = args.BothSeq

    global_Stattext = f"{name}_stat.csv"
    logging.info(f"File prefix name: {name}, Stat file name: {global_Stattext}")
    logging.info("Now, Start transform ...")

    fq1 = f"{name}.FQ1.fastq"
    fq2 = f"{name}.FQ2.fastq"

    command = args.command
    BC = args.BC

    # Define module functions mapping
    module_functions = {
        'all': run_all,
        'transform': run_transform,
        'full': run_full,
        'full2': run_full2,
        'macaca': run_macaca,
        'macaca-transform': run_macaca_transform,
        'simple': run_simple,
        'draw': run_draw,
        'test': run_test,
        'split-migec': run_split_migec,
        'bc-migec-draw': run_bc_migec_draw,
        'migec-draw': run_migec_draw,
        'run3-draw': run_run3_draw,
        'migecfast': run_migecfast
    }

    # Execute the selected module
    if command in module_functions:
        module_functions[command](args, cfg, name, num, TCRtypes, fq1, fq2, global_Stattext)
    else:
        logging.error(f"Unknown command: {command}")
        sys.exit(1)

    # Handle Statindex options if provided
    if args.Statindex:
        handle_statindex(args, cfg, name, num, fq1, global_Stattext)

    endtime = datetime.datetime.now()
    logging.info(f"Total time: {endtime - starttime}")
    logging.info("END")

if __name__ == '__main__':
    main()
