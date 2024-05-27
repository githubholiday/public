file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))/../
ifeq ($(strip $(config)),)
	Bconfig=$(BIN)/config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)
ScriptDir=$(BIN)/bin

Help:
	@echo Description:
	@echo -e "\t" 对T7下机的ATAC数据进行截取和处理，并进行gzip压缩
	@echo Usage:
	@echo -e "\n" Usage获取Config文件
	@echo -e "\t" make -f ${file} project_id= hcsv= outfile= I5_discard= R2_discard=GetSampleList
	@echo  Parameters:
	@echo -e "\t" project_id: 子项目编号
	@echo -e "\t" hcsv: 拆分目录下的Hcsv文件
	@echo -e "\t" I5_discard: 获取I5需要截取掉的长度，使用参数为seqkit -b
	@echo -e "\t" R2_discard: 获取R2需要截取掉的长度，使用参数为seqkit -e
	@echo -e "\t" outfile: 生成的config.ini文件
	@echo -e "\n" Usage:获取ATAC的I5和R2数据长度
	@echo -e "\t" make -f ${file} R2= discard_len= outfile= ATAC_I5
	@echo -e "\t" make -f ${file} R2= discard_len= outfile= ATAC_R2
	@echo  Parameters:
	@echo -e "\t" R2: 样本的R2.fq.gz路径
	@echo -e "\t" discard_len:需要切除的碱基数，I5使用seqkit -b 参数，R2使用seqkit -e参数
	@echo -e "\t" outfile:输出文件名，如outdir/sample1_I5.fq

outdir=$(dir $(abspath $(firstword $(outfile))))
GetSampleList:
	echo "########## GetSampleList start at" `date`
	mkdir -p $(outdir)/
	echo "[Sample]\n" >  $(outfile)
	grep $(project_id) $(hcsv) |awk -F "," '{print $$3"\t"$$11"\t"$(I5_discard)"\t"$(R2_discard)}' >> $(outfile)
	echo "########## GetSampleList  finish at" `date`

outdir=$(dir $(abspath $(firstword $(outfile))))
outfile=$(outdir)/$(sample)/$(sample)_i5.fq
.PHONY:ATAC_Cut
ATAC_I5:
	echo "########## ATAC Get I5 start at" `date`
	mkdir -p $(outdir)/
	$(SEQKIT) trimfq -b $(discard_len) $(R2) > $(outfile)
	$(GZIP) $(outfile)
	echo "########## ATAC Get I5  finish at" `date`

outdir=$(dir $(abspath $(firstword $(outfile))))
outfile=$(outdir)/$(sample)/$(sample)_R2.fq
ATAC_R2:
	echo "########## ATAC Cut R2 start at" `date`
	mkdir $(outdir)
	$(SEQKIT) trimfq -e $(right_bp) $(R2) > $(outfile)
	$(GZIP) $(outfile)
	echo "########## ATAC Cut R2 finish at" `date`