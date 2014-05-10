DATA_CLEAN_DIR=~/repo/thesis/code/tmle
DATA_DUMP_DIR=${DATA_CLEAN_DIR}/data_dump
DB_QUERY=psql --set ON_ERROR_STOP=on -q -d IrisQuebec -f
DB_QUERY_COMMAND=psql --set ON_ERROR_STOP=on -q -d IrisQuebec -c
RSCRIPT=/usr/bin/Rscript --vanilla

all : tables.object

tables.object : ${DATA_DUMP_DIR}/heart.failure.uncombined.object ${DATA_DUMP_DIR}/ami.combined.object ${DATA_DUMP_DIR}/pneumonia.combined.object ${DATA_CLEAN_DIR}/third.R
${RSCRIPT} ${DATA_CLEAN_DIR}/third.R 

${DATA_DUMP_DIR}/heart.failure.combined.object : ${DATA_DUMP_DIR}/heart.failure.uncombined.object ${DATA_CLEAN_DIR}/second.R
${RSCRIPT} ${DATA_CLEAN_DIR}/second.R heart.failure.object

${DATA_DUMP_DIR}/ami.combined.object : ${DATA_DUMP_DIR}/ami.uncombined.object ${DATA_CLEAN_DIR}/second.R
${RSCRIPT} ${DATA_CLEAN_DIR}/second.R ami.object

${DATA_DUMP_DIR}/pneumonia.combined.object : ${DATA_DUMP_DIR}/pneumonia.uncombined.object ${DATA_CLEAN_DIR}/second.R
${RSCRIPT} ${DATA_CLEAN_DIR}/second.R pneumonia.object

${DATA_DUMP_DIR}/heart.failure.uncombined.object : ${DATA_DUMP_DIR}/heart.failure.object ${DATA_CLEAN_DIR}/first.R
${RSCRIPT} ${DATA_CLEAN_DIR}/first.R heart.failure.object

${DATA_DUMP_DIR}/ami.uncombined.object : ${DATA_DUMP_DIR}/ami.object ${DATA_CLEAN_DIR}/first.R
${RSCRIPT} ${DATA_CLEAN_DIR}/first.R ami.object

${DATA_DUMP_DIR}/pneumonia.uncombined.object : ${DATA_DUMP_DIR}/pneumonia.object ${DATA_CLEAN_DIR}/first.R
${RSCRIPT} ${DATA_CLEAN_DIR}/first.R pneumonia.object

# Should be cut up into the three diagnostic categories.
${DATA_DUMP_DIR}/heart.failure.object : ${DATA_DUMP_DIR}/big.matrix.object ${DATA_CLEAN_DIR}/cut.R
${RSCRIPT} ${DATA_CLEAN_DIR}/cut.R heart.failure

# Although it actually outputs several things.
${DATA_DUMP_DIR}/big.matrix.object : ${DATA_DUMP_DIR}/drugs.ahfs.object ${DATA_CLEAN_DIR}/big_matrix.R
${RSCRIPT} ${DATA_CLEAN_DIR}/big_matrix.R

# Build the data files.
${DATA_DUMP_DIR}/drugs.ahfs.object : ${DATA_DUMP_DIR}/quick_tables ${DATA_CLEAN_DIR}/build_data.R
${RSCRIPT} ${DATA_CLEAN_DIR}/build_data.R

# Build the chandan_* prefix tables from aman.readmissions_top20.
# Because SQL doesn't output a file, depend on a empty file in STAMP_DIR
${DATA_DUMP_DIR}/quick_tables : ${DATA_CLEAN_DIR}/quick_tables.sql 
${DB_QUERY} $<
  touch $@
