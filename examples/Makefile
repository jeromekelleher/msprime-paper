# Makefile for the examples.
#

TIME=/usr/bin/time -v

NUM_GWAS_CASES=50000

FILEDIR=files__NOBACKUP__
HDF5_TREES=${FILEDIR}/trees.hdf5
HDF5_TREES_COMPRESSED=${FILEDIR}/trees_compressed.hdf5
TEXT_RECORDS=${FILEDIR}/records.txt
TEXT_RECORDS_COMPRESSED=${FILEDIR}/records.txt.gz
HDF5_MUTATIONS=${FILEDIR}/mutations.hdf5
HDF5_MUTATIONS_COMPRESSED=${FILEDIR}/mutations_compressed.hdf5
NEWICK_SUBSET=${FILEDIR}/newick_subset.txt
NEWICK_SUBSET_COMPRESSED=${FILEDIR}/newick_subset.txt.gz
NEWICK_SUBSET_TREEZIP=${FILEDIR}/newick_subset.txt.trz
NEWICK_SUBSET_TREEZIP_COMPRESSED=${FILEDIR}/newick_subset.txt.trz.gz
SINGLE_NEWICK=${FILEDIR}/single_newick.txt
TREE_BENCHMARKS=${FILEDIR}/tree_benchmarks.txt
TEXT_HAPLOTYPES=${FILEDIR}/haplotypes.txt
TEXT_HAPLOTYPES_COMPRESSED=${FILEDIR}/haplotypes.txt.gz
MACS_OUTPUT=${FILEDIR}/macs.txt
PBWT_HAPLOTYPES=${FILEDIR}/haplotypes.pbwt
PLINK_FILE_PREFIX=${FILEDIR}/plink_input
PLINK_BENCHMARK=${FILEDIR}/plink_benchmarks.txt

HTSLIB=htslib
PBWT=pbwt/pbwt
TREEZIP=treezip/treezip
PLINK_107=p-link
PLINK_19=./plink

all: ${HDF5_TREES_COMPRESSED} ${TEXT_RECORDS_COMPRESSED} \
	${HDF5_MUTATIONS_COMPRESSED} ${NEWICK_SUBSET_COMPRESSED}\
	${TREE_BENCHMARKS} ${TEXT_HAPLOTYPES_COMPRESSED}\
	${PBWT_HAPLOTYPES} ${NEWICK_SUBSET_TREEZIP_COMPRESSED} \
	${PLINK_BENCHMARK}
	
${FILEDIR}:
	mkdir -p ${FILEDIR}

${HDF5}: ${FILEDIR}	
	${TIME} msp simulate 1e5 ${HDF5_TREES} -m 1e8 -r 1e-3 -s 1

${HDF5_TREES_COMPRESSED}: ${HDF5_TREES}
	python examples.py compress ${HDF5_TREES} ${HDF5_TREES_COMPRESSED}

${TEXT_RECORDS}: ${HDF5_TREES}
	msp records ${HDF5_TREES} > ${TEXT_RECORDS}

${TEXT_RECORDS_COMPRESSED}: ${TEXT_RECORDS}
	gzip -c ${TEXT_RECORDS} > ${TEXT_RECORDS_COMPRESSED}

${HDF5_MUTATIONS}: ${HDF5_TREES}
	python examples.py mutate -u 1e-3 ${HDF5_TREES} ${HDF5_MUTATIONS}

${HDF5_MUTATIONS_COMPRESSED}: ${HDF5_MUTATIONS}
	python examples.py compress ${HDF5_MUTATIONS} ${HDF5_MUTATIONS_COMPRESSED}

${NEWICK_SUBSET}: ${HDF5_TREES}
	msp newick ${HDF5_TREES} -p 10 | head -n 1000 > ${NEWICK_SUBSET}

${NEWICK_SUBSET_COMPRESSED}: ${NEWICK_SUBSET}
	gzip -c ${NEWICK_SUBSET} > ${NEWICK_SUBSET_COMPRESSED}

${NEWICK_SUBSET_TREEZIP}: ${NEWICK_SUBSET} ${TREEZIP}
	${TIME} ${TREEZIP} ${NEWICK_SUBSET}

${NEWICK_SUBSET_TREEZIP_COMPRESSED}: ${NEWICK_SUBSET_TREEZIP}
	gzip -c ${NEWICK_SUBSET_TREEZIP} > ${NEWICK_SUBSET_TREEZIP_COMPRESSED}

${SINGLE_NEWICK}: ${NEWICK_SUBSET}
	head -n 1 ${NEWICK_SUBSET} > ${SINGLE_NEWICK}

${TREE_BENCHMARKS}: ${SINGLE_NEWICK}
	python examples.py benchmark-trees ${HDF5_TREES} > ${TREE_BENCHMARKS}
	echo "============\nR Benchmark\n============" >> ${TREE_BENCHMARKS}
	Rscript trees.R >> ${TREE_BENCHMARKS}

${TEXT_HAPLOTYPES}: ${HDF5_MUTATIONS}
	${TIME} msp haplotypes ${HDF5_MUTATIONS} > ${TEXT_HAPLOTYPES}

${TEXT_HAPLOTYPES_COMPRESSED}: ${TEXT_HAPLOTYPES}
	${TIME} gzip -c ${TEXT_HAPLOTYPES} > ${TEXT_HAPLOTYPES_COMPRESSED}

${MACS_OUTPUT}: ${HDF5_MUTATIONS}
	msp macs ${HDF5_MUTATIONS} > ${MACS_OUTPUT}

${PBWT_HAPLOTYPES}: ${MACS_OUTPUT} ${PBWT}
	${TIME} ${PBWT} -readMacs ${MACS_OUTPUT} -write ${PBWT_HAPLOTYPES}

${PLINK_FILE_PREFIX}.ped: ${HDF5_MUTATIONS}
	python examples.py write-ped ${HDF5_MUTATIONS} \
		${PLINK_FILE_PREFIX} --num-cases=${NUM_GWAS_CASES}

${PLINK_FILE_PREFIX}.bed: ${PLINK_FILE_PREFIX}.ped ${PLINK_19}
	${PLINK_19} --file ${PLINK_FILE_PREFIX} --make-bed --out ${PLINK_FILE_PREFIX}

${PLINK_BENCHMARK}: ${PLINK_FILE_PREFIX}.bed ${HDF5_MUTATIONS}
	echo Python > ${PLINK_BENCHMARK}
	${TIME} python examples.py gwas ${HDF5_MUTATIONS} ${FILEDIR}/python-gwas.txt \
		--num-cases=${NUM_GWAS_CASES} 2>> ${PLINK_BENCHMARK}
	echo **plink 1.9** >> ${PLINK_BENCHMARK}
	${TIME} ${PLINK_19} --bfile ${PLINK_FILE_PREFIX} --assoc \
		--out ${FILEDIR}/plink_19 2>> ${PLINK_BENCHMARK}
	echo **plink 1.07** >> ${PLINK_BENCHMARK}
	${TIME} ${PLINK_107} --bfile ${PLINK_FILE_PREFIX} --assoc \
		--out ${FILEDIR}/plink_107 2>> ${PLINK_BENCHMARK}

${HTSLIB}:
	wget https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
	tar -jxf htslib-1.2.1.tar.bz2
	mv htslib-1.2.1 htslib
	${MAKE} -C ${HTSLIB}

${PBWT}: ${HSTLIB}
	git clone https://github.com/richarddurbin/pbwt 
	${MAKE} -C pbwt version.h
	${MAKE} -C pbwt 

${TREEZIP}:
	git clone https://github.com/suzannejmatthews/treezip.git
	cd treezip && ./configure
	make -C treezip
	cp treezip/TreeZip/treezip treezip

${PLINK_19}:
	wget https://www.cog-genomics.org/static/bin/plink151126/plink_linux_x86_64.zip
	unzip plink_linux_x86_64.zip
	rm toy.* prettify LICENSE

clean:
	rm -fR htslib* treezip ${FILEDIR}
