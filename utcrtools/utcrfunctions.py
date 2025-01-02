# utcrtools/functions.py

import os
import re
import gzip
import pandas as pd

def transformer(FRdata, FRdata2, WRdata1, WRdata2, same_seq="CAGTGGTATCAACGCAGAG"):
    while True:
        line1 = FRdata.readline().rstrip().decode('utf-8')
        if not line1:
            break
        line2 = FRdata.readline().rstrip().decode('utf-8')
        lineUMI_before = line2.split(same_seq)
        line3 = FRdata.readline().rstrip().decode('utf-8')
        line4 = FRdata.readline().rstrip().decode('utf-8')
        line2_1 = FRdata2.readline().rstrip().decode('utf-8')
        line2_2 = FRdata2.readline().rstrip().decode('utf-8')
        line2_3 = FRdata2.readline().rstrip().decode('utf-8')
        line2_4 = FRdata2.readline().rstrip().decode('utf-8')

        if len(lineUMI_before) == 2 and len(lineUMI_before[1][0:23]) > 1:
            umi = lineUMI_before[1][0:23]
            if len(umi) == 23 and umi[0] == "T" and umi[5] == "T":
                line4_start = len(lineUMI_before[0])
                line4_end = line4_start + len(umi)
                score = line4[line4_start:line4_end]
                data1 = (line1 + "\tBC UMI:" + umi + ":" + score + "\n" +
                         line2 + "\n" +
                         line3 + "\n" +
                         line4)
                WRdata1.write(data1 + "\n")

                data2 = (line2_1 + "\tSE UMI:" + umi + ":" + score + "\n" +
                         line2_2 + "\n" +
                         line2_3 + "\n" +
                         line2_4)
                WRdata2.write(data2 + "\n")

    WRdata1.close()
    WRdata2.close()

def trans_data(args, name):
    """
    从 args.FQ1 和 args.FQ2 读入 gzip FASTQ，
    找到包含 Sameseq 的 reads 并写到 name.FQ1.fastq / name.FQ2.fastq
    然后再用 FQ2 与 FQ1 交换调用一次
    """
    WR1 = open(name + ".FQ1.fastq", "w+")
    WR2 = open(name + ".FQ2.fastq", "w+")
    transformer(FRdata=gzip.open(args.FQ1, "rb"),
                FRdata2=gzip.open(args.FQ2, "rb"),
                WRdata1=WR1, WRdata2=WR2,
                same_seq=args.BothSeq)

    WR3 = open(name + ".FQ1.fastq", "a")
    WR4 = open(name + ".FQ2.fastq", "a")
    transformer(FRdata=gzip.open(args.FQ2, "rb"),
                FRdata2=gzip.open(args.FQ1, "rb"),
                WRdata1=WR3, WRdata2=WR4,
                same_seq=args.BothSeq)

def migecrun(name, num="3", TCRtypes="TRA,TRB", migecpath=""):
    os.system(f"java -jar {migecpath} Assemble -m {num} -q 15 --mask 0:1 --filter-collisions {name}.FQ1.fastq {name}.FQ2.fastq {name}.cdrdata")
    os.system(f"java -jar {migecpath} CdrBlast -a -R {TCRtypes} {name}.cdrdata/{name}.FQ2.t{num}.cf.fastq {name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt")

def migecrun1(name, num="3", TCRtypes="TRA,TRB", migecpath=""):
    os.system(f"java -jar {migecpath} Assemble -m {num} -q 15 --mask 0:1 --filter-collisions {name}.FQ1.fastq {name}.FQ2.fastq {name}.cdrdata")

def migecrun2(name, num="3", TCRtypes="TRA,TRB", migecpath=""):
    os.system(f"java -jar {migecpath} CdrBlast -a -R {TCRtypes} {name}.cdrdata/{name}.FQ2.t{num}.cf.fastq {name}.cdrdata/{name}_t{num}TRA,TRB.cdrblast.txt")

def migecrun3(name, num="3", TCRtypes="TRA,TRB", SPECIES="HomoSapiens", migecpath=""):
    os.system(f"java -jar {migecpath} CdrBlast -a -R {TCRtypes} {name}.cdrdata/{name}.FQ2.t{num}.cf.fastq {name}.cdrdata/{SPECIES}_{name}_t{num}TRA,TRB.cdrblast.txt")

def migecrunfast(name, num="3", TCRtypes="TRA,TRB", migecpath=""):
    with open("run_migec_cpu.txt", "a+") as migeccpu:
        cmd = f"java -jar {migecpath} Assemble -m {num} -q 15 --mask 0:1 --filter-collisions {name}.FQ1.fastq {name}.FQ2.fastq {name}.cdrdata"
        migeccpu.write(cmd + "\n")

def migecfast(BC2, n=4, TCRtypes="TRA,TRB", migecpath="", PXTCR03_BCsplitpath=""):
    print("Now fast migec")
    os.system("rm -rf run_migec_cpu.txt")

    with open(BC2, "r") as arr:
        for line in arr:
            fqname3 = line.rstrip().split("\t")[0]
            for i in range(1, n):
                num2 = str(i)
                print(i)
                migecrunfast(name=fqname3, num=num2, TCRtypes=TCRtypes, migecpath=migecpath)
    os.system("ParaFly -c run_migec_cpu.txt -CPU 24")

def testfordata(R, endnum=10, Tseq="T....T....T....TCTTGGG"):
    """
    示例函数
    """
    with open(R, "rb") as FR:
        index = 0
        while True:
            line1 = FR.readline().rstrip().decode('utf-8')
            if index == endnum or not line1:
                break
            line2 = FR.readline().rstrip().decode('utf-8')
            line3 = FR.readline().rstrip().decode('utf-8')
            line4 = FR.readline().rstrip().decode('utf-8')
            lineUMI_before = re.findall(Tseq, line2)
            if len(lineUMI_before) != 0:
                print(lineUMI_before[0])
                line2_dist = line2.split(lineUMI_before[0])
                line4_start = len(line2_dist[0])
                line4_end = line4_start + len(lineUMI_before[0])
                score = line4[line4_start:line4_end]
                print(line4_start, line4_end, score)
            index += 1
    return index

def BC_split(BC, fq1, fq2, PXTCR03_BCsplitpath):
    """
    Split BC and run BC_split script
    """
    os.system("rm -rf run_data2.txt run_data2.txt.completed")

    with open(BC, "r") as arr:
        for line in arr:
            fqname = line.rstrip().split("\t")[0] + ".name"
            fqseq = line.rstrip().split("\t")[1]
            print(fqname, fqseq)
            os.system(f"grep -B1 {fqseq} {fq1} | grep @ | awk '{{print \$1}}' > {fqname}")

            # Append commands to run_data2.txt
            os.system(f"echo \"python {PXTCR03_BCsplitpath} {fqname} {fq1}\" >> run_data2.txt")
            os.system(f"echo \"python {PXTCR03_BCsplitpath} {fqname} {fq2}\" >> run_data2.txt")

    os.system("ParaFly -c run_data2.txt -CPU 24")

def count_line(stat_file_name):
    _stat_file_row = os.popen(f"wc -l < {stat_file_name}").read().strip()
    _Fq_data_count = round(int(_stat_file_row) / 4)
    print(f'the reads of {stat_file_name} is {_Fq_data_count}')
    return [_Fq_data_count]

def stat_module(FQone, fq1, stat_file="stat_file", stat_what="p"):
    # fq1 is needed as external variable
    fq_count = count_line(FQone)[0]
    fastq_count = count_line(fq1)[0]
    fastq_pop = f"{fastq_count / fq_count * 100:.2f}%"

    with open(stat_file, "a+") as sf:
        sf.write("File,Reads_count,All_Percentage\n")
        sf.write(f"{FQone},{fq_count},100%\n")
        sf.write(f"{fq1},{fastq_count},{fastq_pop}\n")

    return [fq_count, fastq_count]
