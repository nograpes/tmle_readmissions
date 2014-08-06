CUR_DIR=.
DATA_CLEAN_DIR=${CUR_DIR}/../data_clean
DISEASE_SUBSET_DIR=${CUR_DIR}/disease_subsets

DATA_DUMP_DIR=${CUR_DIR}/data_dump
MATRIX_CACHE_DIR=${CUR_DIR}/matrix_cache

SURVIVAL_DIR=${CUR_DIR}/survival
SURVIVAL_DATA_DUMP_DIR=${SURVIVAL_DIR}/data_dump

TABLES_DIR=${CUR_DIR}/tables
FIGURES_DIR=${CUR_DIR}/figures

RSCRIPT=/usr/bin/Rscript --vanilla
DISEASES=$(wildcard $(DISEASE_SUBSET_DIR)/*)
DISEASES_NO_PATH=$(subst $(DISEASE_SUBSET_DIR)/,,$(DISEASES))

# Without this line, all intermediate files are deleted after make finishes.
# Unfortunately, all targets of implicit rules are intermediate.
# So after the Q model is generated, if Q* fails, then the Q model is deleted.
# Which is unfortunate because it takes two hours for Q to be generated.
.SECONDARY:

all: $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_Q_star_model_%.object)

${FIGURES_DIR}/variable_importance_by_model_and_class.png : $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_G_calibrated_model_%.object)  $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_Q_calibrated_model_%.object) $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/disease_%.object)  ${CUR_DIR}/variable_importance_by_model_and_class.R
	${RSCRIPT} ${CUR_DIR}/variable_importance_by_model_and_class.R

# ${FIGURES_DIR}/tte_distribution.png : $(DISEASES_NO_PATH:%=${SURVIVAL_DATA_DUMP_DIR}/Q_star_survival_%.object) ${SURVIVAL_DIR}/tte_distribution.R
# 	${RSCRIPT} ${SURVIVAL_DIR}/tte_distribution.R

# $(DISEASES_NO_PATH:%=${SURVIVAL_DATA_DUMP_DIR}/Q_star_survival_%.object)
${TABLES_DIR}/disease.results.table.object : $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/rf_Q_star_model_%.object)	 $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/disease_%.object)	    $(DISEASES_NO_PATH:%=${DATA_DUMP_DIR}/crude_readmissions_risk_%.object)	
	${RSCRIPT} ${CUR_DIR}/disease_tables.R $@

${DATA_DUMP_DIR}/rf_Q_star_model_%.object :  ${DATA_DUMP_DIR}/rf_G_calibrated_model_%.object ${DATA_DUMP_DIR}/rf_Q_calibrated_model_%.object ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_Q_star_model.R
	${RSCRIPT} ${CUR_DIR}/build_rf_Q_star_model.R $^ $@ ${MATRIX_CACHE_DIR}

# Crude readmission risk
${DATA_DUMP_DIR}/crude_readmissions_risk_%.object	: ${DATA_DUMP_DIR}/disease_%.object  ${CUR_DIR}/crude_readmission_risk.R
		${RSCRIPT} ${CUR_DIR}/crude_readmission_risk.R $< $@
		
# Survival targets.	
# ${SURVIVAL_DATA_DUMP_DIR}/Q_star_survival_ami.object : ${DATA_DUMP_DIR}/disease_ami.object ${DATA_DUMP_DIR}/rf_G_calibrated_model_ami.object ${SURVIVAL_DATA_DUMP_DIR}/glmnet_g_censor_ami.object ${SURVIVAL_DATA_DUMP_DIR}/glmnet_Q_ami.object ${SURVIVAL_DIR}/Q_star_survival.R
# 	${RSCRIPT} ${SURVIVAL_DIR}/Q_star_survival.R $^ $@

${SURVIVAL_DATA_DUMP_DIR}/Q_star_survival_%.object : ${DATA_DUMP_DIR}/disease_%.object ${DATA_DUMP_DIR}/rf_G_calibrated_model_%.object ${SURVIVAL_DATA_DUMP_DIR}/glmnet_g_censor_%.object ${SURVIVAL_DATA_DUMP_DIR}/glmnet_Q_%.object ${SURVIVAL_DIR}/Q_star_survival.R
	${RSCRIPT} ${SURVIVAL_DIR}/Q_star_survival.R $^ $@
	
${SURVIVAL_DATA_DUMP_DIR}/glmnet_g_censor_%.object : ${DATA_DUMP_DIR}/disease_%.object ${SURVIVAL_DIR}/glmnet_g_censor.R
	${RSCRIPT} ${SURVIVAL_DIR}/glmnet_g_censor.R $< $@

${SURVIVAL_DATA_DUMP_DIR}/glmnet_Q_%.object : ${DATA_DUMP_DIR}/disease_%.object ${SURVIVAL_DIR}/glmnet_Q.R
	${RSCRIPT} ${SURVIVAL_DIR}/glmnet_Q.R $< $@
# End of survival targets.
	
#${DATA_DUMP_DIR}/glmnet_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_glmnet_Q_model.R
#	${RSCRIPT} ${CUR_DIR}/build_glmnet_Q_model.R $< $@

# ${DATA_DUMP_DIR}/rf_Q_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_Q_model.R
# 	${RSCRIPT} ${CUR_DIR}/build_rf_Q_model.R $< $@ ${MATRIX_CACHE_DIR}

# ${DATA_DUMP_DIR}/rf_G_calibrated_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_G_calibrated_model.R
# 	${RSCRIPT} ${CUR_DIR}/build_rf_G_calibrated_model.R $< $@ ${MATRIX_CACHE_DIR}

${DATA_DUMP_DIR}/rf_Q_calibrated_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_Q_model.R
	${RSCRIPT} ${CUR_DIR}/build_rf_Q_calibrated_model.R $< $@ ${MATRIX_CACHE_DIR}

${DATA_DUMP_DIR}/rf_G_calibrated_model_%.object : ${DATA_DUMP_DIR}/disease_%.object ${CUR_DIR}/build_rf_G_calibrated_model.R
	${RSCRIPT} ${CUR_DIR}/build_rf_G_calibrated_model.R $< $@ ${MATRIX_CACHE_DIR}
	
# This idiom is the only way I know to have a pattern depend on a shell wildcard
# where you need to pass the shell wildcard as an argument, but not the other dependencies.
${DATA_DUMP_DIR}/disease_%.object : ${DISEASES}
	${RSCRIPT} ${CUR_DIR}/build_data.R ${DATA_CLEAN_DIR}/data_source.R ${DATA_DUMP_DIR} $^

# Other dependencies.
${DISEASES} : ${CUR_DIR}/build_data.R 
	touch $@
# End of idiom

${CUR_DIR}/build_data.R : ${DATA_CLEAN_DIR}/data_source.R ${DATA_DUMP_DIR}/tables_finished_stamp
	
${DATA_DUMP_DIR}/tables_finished_stamp : ${CUR_DIR}/quick_tables.sql
	psql -q -d IrisQuebec -f $<
	touch $@