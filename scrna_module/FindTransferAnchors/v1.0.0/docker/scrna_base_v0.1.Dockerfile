FROM komais/centos_base_full:v0.0.1

#维护者信息
LABEL author="leiguo"
USER root

#Installing Monocle 3
RUN	yum install -y cmake && \
	conda install -y -c conda-forge imagemagick && \
	conda install -y -c conda-forge r-devtools && \
	conda install -y -c conda-forge r-terra && \
	conda install -y -c conda-forge r-ggrastr && \
	conda install -y -c conda-forge r-seurat && \
	conda install -y -c bioconda r-monocle3 && \
	conda install -y -c bioconda r-harmony && \
	conda install -y -c bioconda bioconductor-monocle && \
	conda install -y -c bioconda bioconductor-singler


#scibet
RUN	Rscript -e 'install.packages("tidyverse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' &&\
	Rscript -e 'devtools::install_github("PaulingLiu/scibet")' && \
	Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")' && \
	Rscript -e 'install.packages("ggsignif", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' && \
	Rscript -e 'install.packages("R.utils", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' && \
	Rscript -e 'install.packages("configr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' && \
	Rscript -e 'install.packages("argparser", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")'


RUN	Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.tar.gz", repos=NULL, type="source")' && \
	Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.3.tar.gz", repos=NULL, type="source")' && \
	Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.3.0.1.tar.gz", repos=NULL, type="source")' && \
	Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.5.0.tar.gz", repos=NULL, type="source")' && \
	Rscript -e 'install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/leidenbase_0.1.25.tar.gz", repos=NULL, type="source")'


RUN	Rscript -e 'devtools::install_github("immunogenomics/harmony")'

USER test_user
