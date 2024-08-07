file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
ifdef config
	include $(config)
else
	include $(BIN)/config/config.txt
endif
include $(Bconfig)
ScriptDir=$(BIN)/bin

HELP:
	@echo Description: 此脚本用于ccs reads比对
	@echo Usage :
	@echo make -f makefile config= ref= sample= cleandir= indir= outdir= Build_Index Alignment RM_Host Host_Stat
	@echo
	@echo [参数说明]    :
	@echo config        : [可选]流程配置文件，默认为本目录下config/config.txt
	@echo ref           : [必选]宿主参考基因组.fa文件
	@echo sample        : [必选]样品名称,输出文件前缀
	@echo R1            : [必选]比对R1.fq.gz路径
	@echo R2            : [必选]比对R2.fq.gz路径
	@echo outdir        : [必选]输出目录
	@echo
	@echo [Version]:
	@echo v1.0.0 2021-11-15 by 'Simon Lee' [huayunli\@genome.cn];

BWA_Index:
	echo Index start at `date`
	$(BWA) index $(ref)
	echo Index finished at `date`


Alignment:
	echo Alignment_gene start at `date`
	[ -s ${outdir} ] || mkdir -p ${outdir}
	$(BWA) mem $(ref) $(R1) $(R2) -t 8 | $(SAMTOOLS) view -bS - > ${outdir}/${sample}.bam && echo "done!"
	${SAMTOOLS} sort -@ 4 -m 2G ${outdir}/${sample}.bam -o ${outdir}/${sample}.sort.bam && ${SAMTOOLS} flagstat ${outdir}/${sample}.sort.bam > ${outdir}/${sample}.stat.xls && echo "done!"
	echo Alignment_gene finished at `date`

Host_RM_PE:
	echo Alignment start at `date`
	[ -s ${outdir} ] || mkdir -p ${outdir}
	if [ $(ref) == "no_ref" ] ;\
	then \
		ln -s ${R1} ${outdir}/${sample}_umap_1.fq.gz ;\
		ln -s ${R2} ${outdir}/${sample}_umap_2.fq.gz ;\
	else \
		$(BWA) mem $(ref) ${R1} ${R2} -t 8 | $(SAMTOOLS) view -bS - > ${outdir}/${sample}.aln.bam ;\
		${SAMTOOLS} sort -@ 4 -m 2G ${outdir}/${sample}.aln.bam -o ${outdir}/${sample}.aln.sort.bam ;\
		${SAMTOOLS} flagstat ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.stat.xls ;\
		${SAMTOOLS} view -b -f 4 ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.unmapped.bam ;\
		${BAM2FASTQ} -f -o ${outdir}/${sample}_umap#.fq ${outdir}/${sample}.unmapped.bam ;\
		$(GZIP) ${outdir}/${sample}_umap_1.fq ;\
		$(FQTOOLS) stat ${outdir}/${sample}_umap_1.fq.gz ${outdir}/${sample}_umap_1.fq.gz.report ;\
		$(GZIP) ${outdir}/${sample}_umap_2.fq ;\
		$(FQTOOLS) stat ${outdir}/${sample}_umap_2.fq.gz ${outdir}/${sample}_umap_2.fq.gz.report ;\
		$(PYTHON3) $(BIN)/script/stat_host.py -i ${outdir}/${sample}.stat.xls -r  ${outdir}/${sample}_umap_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls ;\
	fi
	echo Alignment finished at `date`

Host_RM_SE:
	echo Alignment start at `date`
	[ -s ${outdir} ] || mkdir -p ${outdir}
	if [ $(ref) == "no_ref" ] ;\
	then \
		ln -s ${R1} ${outdir}/${sample}_umap_1.fq.gz ;\
		ln -s ${R2} ${outdir}/${sample}_umap_2.fq.gz ;\
	else \
		$(BWA) mem $(ref) ${R1} -t 8 | $(SAMTOOLS) view -bS - > ${outdir}/${sample}.aln.bam ;\
		${SAMTOOLS} sort -@ 4 -m 2G ${outdir}/${sample}.aln.bam -o ${outdir}/${sample}.aln.sort.bam ;\
		${SAMTOOLS} flagstat ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.stat.xls ;\
		${SAMTOOLS} view -b -f 4 ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.unmapped.bam ;\
		${BAM2FASTQ} -f -o ${outdir}/${sample}_umap_1.fq ${outdir}/${sample}.unmapped.bam ;\
		$(GZIP) ${outdir}/${sample}_umap_1.fq ;\
		$(FQTOOLS) stat ${outdir}/${sample}_umap_1.fq.gz ${outdir}/${sample}_umap_1.fq.gz.report ;\
		$(PYTHON3) $(BIN)/script/stat_host.py -i ${outdir}/${sample}.stat.xls -r  ${outdir}/${sample}_umap_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls ;\
	fi
	echo Alignment finished at `date`

BWA_MEM:
	echo bwa mem start at `date`
	$(BWA) mem $(ref) ${R1} ${R2} -t 8 | $(SAMTOOLS) view -bS - > ${outdir}/${sample}.aln.bam
	${SAMTOOLS} sort -@ 4 -m 2G ${outdir}/${sample}.aln.bam -o ${outdir}/${sample}.aln.sort.bam 
	${SAMTOOLS} flagstat ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.stat.xls
	${SAMTOOLS} view -b -F 4 ${outdir}/${sample}.aln.sort.bam > ${outdir}/${sample}.mapped.bam
	echo bwa mem end at `date`

BAM2FQ_SE:
	echo bam2fq SE start at `date`
	${BAM2FASTQ} -f -o ${outdir}/${sample}_umap_1.fq ${outdir}/${sample}.unmapped.bam
	$(GZIP) ${outdir}/${sample}_umap_1.fq
	$(FQTOOLS) stat ${outdir}/${sample}_umap_1.fq.gz ${outdir}/${sample}_umap_1.fq.gz.report
	$(PYTHON3) $(BIN)/script/stat_host.py -i ${outdir}/${sample}.stat.xls -r  ${outdir}/${sample}_umap_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls
	echo bam2fq SE end at `date`

BAM2FQ_PE:
	echo bam2fq PE start at `date`
	${BAM2FASTQ} -f -o ${outdir}/${sample}_umap_#.fq ${outdir}/${sample}.unmapped.bam
	$(GZIP) ${outdir}/${sample}_umap_1.fq
	$(GZIP) ${outdir}/${sample}_umap_2.fq
	$(FQTOOLS) stat ${outdir}/${sample}_umap_1.fq.gz ${outdir}/${sample}_umap_1.fq.gz.report
	$(FQTOOLS) stat ${outdir}/${sample}_umap_2.fq.gz ${outdir}/${sample}_umap_2.fq.gz.report
	$(PYTHON3) $(BIN)/script/stat_host.py -i ${outdir}/${sample}.stat.xls -r  ${outdir}/${sample}_umap_1.fq.gz.report -s ${sample} -o ${outdir}/${sample}.map.stat.xls
	echo bam2fq PE end at `date`

Host_Stat:
	echo Host_Stat start at `date`
	if [ $(ref) == "no_ref" ] ;\
	then \
		echo no ref ;\
	else \
		$(CSVTK) -t join  -i ${infile} -o ${outdir}/Align_stat.xls  ;\
	fi
	echo Host_Stat finished at `date`

QC:
	echo QC start at `date`
	if [ -f $(filter_stat) ]; then\
		${PERL} $(BIN)/script/filter_rename.pl $(filter_stat) ${info} >$(outdir)/STAT_result.xls ;\
	fi
	cat ${info}|awk -F "\t" '{print $$4 "\t" $$3 "\t" $$5}'|sed '1iSample\tDescription\tGroup' >$(outdir)/sample_info.xls
	echo QC finished at `date`

Count:
	echo Stat all sample start at `date`
	$(PYTHON3) $(BIN)/script/gene_count.py -i $(input) -o $(output) -s $(SAMTOOLS)
	echo Stat all sample end at `date`

RM_Count:
	echo RM_Count start at `date`
	mkdir -p $(outdir)
	$(PYTHON3) $(BIN)/script/filter_count.py -i $(input) -o $(output) -n $(gene_count)
	${PERL} $(BIN)/script/fa_get_list.pl $(fasta) $(output) >$(outfasta)
	echo RM_Count end at `date`

