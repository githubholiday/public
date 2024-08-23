file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
ifdef config
	include $(config)
else
	include $(BIN)/config/config.txt
endif
include $(Bconfig)
ScriptDir=$(BIN)/script
FQTOOLS=$(ScriptDir)/fqtools

HELP:
	@echo Description: 此脚本用于ccs reads比对
	@echo Usage :
	@echo make -f makefile config= ref= BWA_Index
	@echo   [参数说明]
	@echo   用于ref建库,ref为fa文件
	@echo make -f makefile config= ref= R1= R2= outdir= sample= BWA_MEM
	@echo   [参数说明]
	@echo   ref:参考基因组fa文件(建库完成的,同级目录下有建库文件)
	@echo   R1/R2:用于比对的R1和R2的fq文件,如果无R2,可以不给
	@echo   outdir:输出目录,也可以使用aln_bam和sort_bam参数
	@echo   sample:输出文件前缀,也可以使用aln_bam和sort_bam参数
	@echo   flag_stat_file:samtools flagstat输出文件
	@echo   输出:outdir/sample.aln.bam和outdir/sample.aln.sort.bam
	@echo make -f makefile config= inbam= outbam= para= SelectBam
	@echo   [参数说明]
	@echo   获取比对上/未比对上的bam
	@echo   inbam:sort后的bam文件
	@echo   outbam:提取后的bam文件
	@echo   para:参数,可选[-F,-f],其中-F为提取比对上的bam,-f为提取未比对上的bam
	@echo make -f makefile config= inbam= flag_stat_file= depth_file= coverage_out= sample= DepthStat
	@echo   [参数说明]	
	@echo   对sort后的bam进行flagstat和深度/覆盖度统计
	@echo   inbam:sort后的bam文件
	@echo   flag_stat_file:输入-flagstat后的输出文件
	@echo   depth_file:genomeCoverageBed的输出文件
	@echo   coverage_out:覆盖度统计输出文件
	@echo   sample:样本名称
	@echo make -f makefile config= inbam= sample= type= flag_stat_file= outdir= BAM2FQ_SE
	@echo   [参数说明]	
	@echo   将输入的bam转化为SE的fq文件,并且统计比对上和未比对上reads数/比率
	@echo   inbam:sort后的bam文件
	@echo   sample:样本名称
	@echo   type:输入的bam是比对上的还是未比对上的，map|unmap
	@echo   flag_stat_file:flagstat后的输出文件
	@echo   outdir:输出目录
	@echo   输出信息：outdir/sample_type_1|2.fq|fq.gz,outdir/sample_type_1.fq.gz.report,outdir/sample_map.stat.xls,
	@echo   sample:样本名称
	@echo make -f makefile config= inbam= sample= type= flag_stat_file= outdir= BAM2FQ_PE
	@echo   [参数说明]	
	@echo   将输入的bam转化为SE的fq文件,并且统计比对上和未比对上reads数/比率
	@echo   inbam:sort后的bam文件
	@echo   sample:样本名称
	@echo   type:输入的bam是比对上的还是未比对上的，map|unmap
	@echo   flag_stat_file:flagstat后的输出文件
	@echo   outdir:输出目录
	@echo   输出信息：outdir/sample_type_1|2.fq|fq.gz,outdir/sample_type_1.fq.gz.report,outdir/sample_map.stat.xls,
	@echo   sample:样本名称

BWA_Index:
	echo Index start at `date`
	$(BWA) index $(ref)
	echo Index finished at `date`

aln_bam=${outdir}/${sample}.aln.bam
sort_bam=${outdir}/${sample}.aln.sort.bam
BWA_MEM:
	echo bwa mem and sort start at `date`
	mkdir -p $(outdir)
	$(BWA) mem $(ref) ${R1} ${R2} -t 8 | $(SAMTOOLS) view -bS - > ${aln_bam}
	${SAMTOOLS} sort -@ 4 -m 2G ${aln_bam} -o ${sort_bam}
	${SAMTOOLS} flagstat ${sort_bam} > ${flag_stat_file}
	echo bwa mem end and sort at `date`

outdir=$(dir $(abspath $(firstword $(outbam))))
para=-F #比对上的，-f未比对上的
SelectBam:
	echo get map or unmap bam from sort.bam start at `date`
	mkdir -p $(outdir)
	${SAMTOOLS} view -b ${para} 4 ${inbam} > ${outbam}
	echo get map or unmap bam from sort.bam end at `date`

DepthStat:
	echo bam stat start at `date`
	genomeCoverageBed -ibam ${inbam} > $(depth_file)
	$(PYTHON3) $(ScriptDir)/coverage.py -f $(flag_stat_file) -d $(depth_file) -o $(coverage_out) -n $(sample)

prefix=@$(sample)_$(type)
BAM2FQ_SE:
	echo bam2fq SE start at `date`
	mkdir -p $(outdir)
	${BAM2FASTQ} -f -o ${outdir}/${prefix}_1.fq ${inbam}
	$(GZIP) ${outdir}/${prefix}_1.fq
	$(FQTOOLS) stat ${outdir}/${prefix}_1.fq.gz ${outdir}/${prefix}_1.fq.gz.report
	$(PYTHON3) $(ScriptDir)/stat_host.py -i ${flag_stat_file} -r ${outdir}/${prefix}_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls -t $(type)
	echo bam2fq SE end at `date`

prefix=$(sample)_$(type)
BAM2FQ_PE:
	echo bam2fq PE start at `date`
	mkdir -p $(outdir)
	${BAM2FASTQ} -f -o ${outdir}/${prefix}_#.fq ${inbam}
	$(GZIP) ${outdir}/${prefix}_1.fq
	$(GZIP) ${outdir}/${prefix}_2.fq
	$(FQTOOLS) stat ${outdir}/${prefix}_1.fq.gz ${outdir}/${prefix}_1.fq.gz.report
	$(FQTOOLS) stat ${outdir}/${prefix}_2.fq.gz ${outdir}/${prefix}_2.fq.gz.report
	$(PYTHON3) $(ScriptDir)/script/stat_host.py -i ${flag_stat_file} -r  ${outdir}/${prefix}_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls -t $(type)
	echo bam2fq PE end at `date`

