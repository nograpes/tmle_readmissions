DATA_CLEAN_DIR=.
DATA_DUMP_DIR=${DATA_CLEAN_DIR}/data_dump
DISEASE_SUBSET_DIR=${DATA_CLEAN_DIR}/disease_subsets
RSCRIPT=/usr/bin/Rscript --vanilla
DISEASES=$(wildcard $(DISEASE_SUBSET_DIR)/*)
DISEASES_NO_PATH=$(subst $(DISEASE_SUBSET_DIR)/,,$(DISEASES))

# Without this line, all intermediate files are deleted after make finishes.
# Unfortunately, all targets of implicit rules are intermediate.
# So after the Q model is generated, if Q* fails, then the Q model is deleted.
# Which is unfortunate because it takes two hours for Q to be generated.
.SECONDARY:

all: $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_Q_star_model_%.object)

${DATA_DUMP_DIR}/rf_Q_star_model_%.object : ${DATA_DUMP_DIR}/rf_G_model_%.object ${DATA_DUMP_DIR}/rf_Q_model_%.object ${DATA_DUMP_DIR}/disease_%.object  ${DATA_CLEAN_DIR}/build_rf_Q_star_model.R 
	${RSCRIPT} ${DATA_CLEAN_DIR}/build_rf_Q_star_model.R ${DATA_DUMP_DIR}/rf_G_model_$*.object ${DATA_DUMP_DIR}/rf_Q_model_$*.object ${DATA_DUMP_DIR}/disease_$*.object $@
	
${DATA_DUMP_DIR}/rf_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${DATA_CLEAN_DIR}/build_rf_Q_model.R
	${RSCRIPT} ${DATA_CLEAN_DIR}/build_rf_Q_model.R $< $@

${DATA_DUMP_DIR}/rf_G_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${DATA_CLEAN_DIR}/build_rf_G_model.R
	${RSCRIPT} ${DATA_CLEAN_DIR}/build_rf_G_model.R $< $@

# This idiom is the only way I know to have a pattern depend on a shell wildcard
# where you need to pass the shell wildcard as an argument, but not the other dependencies.
${DATA_DUMP_DIR}/disease_%.object : ${DISEASES}
	${RSCRIPT} ${DATA_CLEAN_DIR}/cut.R $^

# Other dependencies.
${DISEASE_SUBSET_DIR}/* : ${DATA_CLEAN_DIR}/cut.R ${DATA_DUMP_DIR}/fixed.vars.mat.object
	touch $@
# End of idiom

${DATA_DUMP_DIR}/fixed.vars.mat.object : ${DATA_CLEAN_DIR}/build_data.R
	${RSCRIPT} $< $@

