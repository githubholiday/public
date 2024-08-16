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
	@echo -e "\t" make -f ${file} R1= outdir= sample= cpu= other_parameter= FastpSE
	@echo -e "\t" make -f ${file} R1= R2= outdir= sample= cpu= other_parameter= FastpPE
	@echo  Parameters:
	@echo -e "\t" R1/R2: 过滤输入的R1和R2的fq
	@echo -e "\t" sample:输出的样本名称
	@echo -e "\t" outdir:输出目录
	@echo -e "\t" cpu:使用的CPU数量，默认为3
	@echo -e "\t" make -f ${file} input= output= FilterStat
	@echo  Parameters:
	@echo -e "\t" input: 过滤后的json文件，可以给成 */*.json文件
	@echo -e "\t" output：整理输出文件

FastpSE:
	echo filter SE start at `date`
	mkdir -p $(outdir)
	$(fastp) --thread ${cpu} -c --in1 ${R1} --out1 ${outdir}/${sample}_R1.fq.gz -j ${outdir}/${sample}.stat.json -h ${outdir}/${sample}.stat.html ${other_parameter} 
	echo filter SE finished at `date`

FastpPE:
	echo filter SE start at `date`
	mkdir -p $(outdir)
	$(fastp) --thread ${cpu} -c --in1 ${R1} --out1 ${outdir}/${sample}_R1.fq.gz --in2 ${R2} --out2 ${outdir}/${sample}_R2.fq.gz -j ${outdir}/${sample}.stat.json -h ${outdir}/${sample}.stat.html ${other_parameter} 
	echo filter SE finished at `date`

outdir=$(dir $(abspath $(firstword $(output))))
.PHONY:FilterStat
FilterStat:
	echo "########## filter_stat  start at" `date`
	mkdir -p $(outdir)
	$(BIN)/merge_json --input $(input) --output $(output)
	echo "########## filter_stat finish at" `date`