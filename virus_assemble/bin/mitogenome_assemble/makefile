file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
ifdef config
	include $(config)
else
	include $(BIN)/config/config.txt
endif
include $(Bconfig)
ScriptDir=$(BIN)/script

cpu=?3
fastp=fastp
HELP:
	@echo Description: 此脚本用于ccs reads比对
	@echo Usage :
	@echo -e "\t" make -f ${file} R1= R2= genome_type= max_round= outdir=  GetOrangellePE
	@echo -e "\t" make -f ${file} R1= genome_type= max_round= outdir= GetOrangelleSE
	@echo  Parameters:
	@echo -e "\t" R1/R2: 过滤输入的R1和R2的fq
	@echo -e "\t" outdir:输出目录
	@echo -e "\t" genome_type:使用的CPU数量，默认为3
	@echo -e "\t" max_round:使用的CPU数量，默认为3
	@echo -e "\t" make -f ${file} input= output= FilterStat
	@echo  Parameters:
	@echo -e "\t" input: 过滤后的json文件，可以给成 */*.json文件
	@echo -e "\t" output：整理输出文件

GetOrangellePE:
	echo GetOrangellePE start at `date`
	mkdir -p $(outdir)
	get_organelle_from_reads.py -1 $(R1) -2 $(R2) -F $(genome_type) -R $(max_round) -o $(outdir) 
	echo GetOrangellePE finished at `date`

GetOrangelleSE:
	echo GetOrangelleSE start at `date`
	mkdir -p $(outdir)
	get_organelle_from_reads.py -u $(R1) -F $(genome_type) -R $(max_round) -o $(outdir) 
	echo GetOrangelleSE finished at `date`

Accesment:
	echo Accesment start at `date`
	mkdir -p $(outdir)
	grep circular -A 1 $(infa) > $(out_fa)
	[ -d $(out_fa) ] && quast.py $(out_fa) -r $(ref) || echo no circular fa
	echo Accesment finished at `date`	