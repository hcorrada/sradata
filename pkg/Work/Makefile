RUNIT_DIR=$(shell pwd)
PKG_PATH=$(shell dirname ${RUNIT_DIR})

SOURCES=$(wildcard ${PKG_PATH}/src/*.c)
RSOURCES=$(wildcard ${PKG_PATH}/R/*.R)

SRALIB=${HOME}/sgseq/src/intmodel/sra/lib64
SRAILIB=${HOME}/sgseq/src/intmodel/sra/ilib64
SRA_PATH=${HOME}/sgseq/src/intmodel/sra/sra_sdk-1.0.0-b10

TBLPATH=/amber1/archive/sgseq/data/1kgenomes/low_coverage/HG00116
ACCESSION=SRR035464


all: install test
test: install
	Rscript --slave test.R ${TBLPATH}/${ACCESSION}

${PKG_PATH}/configure: ${PKG_PATH}/configure.ac
	cd ${PKG_PATH}; \
	autoconf

install: ${PKG_PATH}/configure ${SOURCES} ${RSOURCES}
	R CMD INSTALL --configure-args="--with-sra-includes=${SRA_PATH}/itf --with-sra-libs=${SRALIB}" ${PKG_PATH}

installSite: ${PKG_PATH}/configure ${SOURCES} ${RSOURCES}
	R --args '--with-sra-includes=${SRA_PATH}/itf --with-sra-libs=${SRALIB}'

.PHONY: all install test
