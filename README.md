(py38)utcrtools -h
usage: utcrtools [-h] --FQ1 FQ1 --FQ2 FQ2 --BC BC --ini INI [--FN FN] [--TCR TCR] [--Num NUM] [--BothSeq BOTHSEQ] [--Statwhat {a,a1,p}] [--Statindex STATINDEX] [--DF DF] [--SV SV]
                 {all,transform,full,full2,macaca,macaca-transform,simple,draw,test,split-migec,bc-migec-draw,migec-draw,run3-draw,migecfast} ...

UTCRtools: TCR processing suite.

positional arguments:
  {all,transform,full,full2,macaca,macaca-transform,simple,draw,test,split-migec,bc-migec-draw,migec-draw,run3-draw,migecfast}
                        Available modules
    all                 Transform + split BC + run Migec three times
    transform           Only transform data
    full                Transform + Migec + draw
    full2               Transform + run Migec + draw
    macaca              Run Migec for MacacaMulatta
    macaca-transform    Transform + run Migec for MacacaMulatta
    simple              Simple Migec run + draw
    draw                Only draw
    test                Test data
    split-migec         Transform + BC_split + Migec
    bc-migec-draw       BC_split + Migec + draw
    migec-draw          Directly run Migec + draw on split files
    run3-draw           Run Migec three times on split files + draw
    migecfast           Migecfast + run Migec2 + draw

optional arguments:
  -h, --help            show this help message and exit
  --FQ1 FQ1             Fastq R1
  --FQ2 FQ2             Fastq R2
  --BC BC               File with BC name and BC seq
  --ini INI             config.ini
  --FN FN               Filename prefix
  --TCR TCR             TCR types:<TRA|TRB|TRA,TRB|...>, default:TRA,TRB
  --Num NUM             TCR num, default:3
  --BothSeq BOTHSEQ     Same sequence in your data, default:CAGTGGTATCAACGCAGAG
  --Statwhat {a,a1,p}   Data to stat: a (all), a1 (some), p (primary data)
  --Statindex STATINDEX
                        Plot and stat folder. Run Statwhat module before if using this.
  --DF DF               Folder to get statistics data
  --SV SV               Folder to save statistics data
