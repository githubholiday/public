FROM komais/centos_base_full:v0.0.1

#维护者信息
LABEL author="leiguo"
USER root

#软件绝对路径
#R:/software/conda/bin/R
#python3:/software/conda/bin/python3
#convert:/software/conda/bin/convert


#Installing spotlight
RUN	conda install -y -c conda-forge imagemagick && \
	conda install -y -c conda-forge r-seurat r-argparse r-devtools && \
	conda install -y -c bioconda bioconductor-spotlight && \
	conda install -y -c bioconda bioconductor-scran bioconductor-scater && \
	conda install -y -c bioconda bioconductor-singlecellexperiment


#scibet
RUN	Rscript -e 'install.packages("configr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' && \
	Rscript -e 'install.packages("ggnewscale", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' && \
	Rscript -e 'install.packages("tidyverse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")'


RUN	Rscript -e 'devtools::install_github("PaulingLiu/scibet")'


USER test_user