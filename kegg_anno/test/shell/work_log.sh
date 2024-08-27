#第一步：整理kegg结果，并将tpm结果追加到注释结果后面
make -f ../../mk_kegg input=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/input/*.out kegg_level_file=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/KEGG/current/data/map_pathway.list TPM=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/input/TPM.all.xls outdir=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result/ KEGG_Combine
#第二步：计算KEGG每个level的表达量信息
make -f ../../mk_kegg infile=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result//Gene_All.KEGG.TPM.xls outdir=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result KEGG_level_tpm_stat
#第三步：整理交付目录
make -f ../../mk_kegg indir=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result/ outdir=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result/ upload

#合并运行 
make -f ../../mk_kegg input=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/input/*.out kegg_level_file=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/KEGG/current/data/map_pathway.list TPM=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/input/TPM.all.xls outdir=/annogene/data1/bioinfo/Seq/RD/PMO/tuchengfang/238/Develop/public/kegg_anno/test/result/ KEGG_Combine KEGG_level_tpm_stat upload
