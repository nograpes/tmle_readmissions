CUR_DIR=.
DATA_CLEAN_DIR=${CUR_DIR}/../data_clean
DISEASE_SUBSET_DIR=${CUR_DIR}/disease_subsets

DATA_DUMP_DIR=${CUR_DIR}/data_dump
MATRIX_CACHE_DIR=${CUR_DIR}/matrix_cache

RSCRIPT=/usr/bin/Rscript --vanilla
DISEASES=$(wildcard $(DISEASE_SUBSET_DIR)/*)
DISEASES_NO_PATH=$(subst $(DISEASE_SUBSET_DIR)/,,$(DISEASES))

# Without this line, all intermediate files are deleted after make finishes.
# Unfortunately, all targets of implicit rules are intermediate.
# So after the Q model is generated, if Q* fails, then the Q model is deleted.
# Which is unfortunate because it takes two hours for Q to be generated.
.SECONDARY:

all: $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_Q_star_model_%.object)

${DATA_DUMP_DIR}/rf_Q_star_model_%.object : ${DATA_DUMP_DIR}/rf_G_model_%.object ${DATA_DUMP_DIR}/rf_Q_model_%.object ${DATA_DUMP_DIR}/glm_Q_model_%.object ${DATA_DUMP_DIR}/glmnet_Q_model_%.object ${DATA_DUMP_DIR}/disease_%.object  ${CUR_DIR}/build_rf_Q_star_model.R 
	${RSCRIPT} ${CUR_DIR}/build_rf_Q_star_model.R ${DATA_DUMP_DIR}/rf_G_model_$*.object ${DATA_DUMP_DIR}/rf_Q_model_$*.object ${DATA_DUMP_DIR}/disease_$*.object $@
	
${DATA_DUMP_DIR}/glm_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_glm_Q_model.R
	${RSCRIPT} ${CUR_DIR}/build_glm_Q_model.R $< $@

${DATA_DUMP_DIR}/glmnet_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_glmnet_Q_model.R
	${RSCRIPT} ${CUR_DIR}/build_glmnet_Q_model.R $< $@

${DATA_DUMP_DIR}/rf_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_Q_model.R
	${RSCRIPT} ${CUR_DIR}/build_rf_Q_model.R $< $@ ${MATRIX_CACHE_DIR}

${DATA_DUMP_DIR}/rf_G_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_G_model.R
	${RSCRIPT} ${CUR_DIR}/build_rf_G_model.R $< $@ ${MATRIX_CACHE_DIR}

# This idiom is the only way I know to have a pattern depend on a shell wildcard
# where you need to pass the shell wildcard as an argument, but not the other dependencies.
${DATA_DUMP_DIR}/disease_%.object : ${DISEASES}
	${RSCRIPT} ${CUR_DIR}/build_data.R ${DATA_CLEAN_DIR}/data_source.R ${DATA_DUMP_DIR} $^

# Other dependencies.
${DISEASES} : ${CUR_DIR}/build_data.R 
	touch $@
# End of idiom

${CUR_DIR}/build_data.R : ${DATA_CLEAN_DIR}/data_source.R
