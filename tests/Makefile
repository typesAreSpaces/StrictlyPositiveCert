QUIET_MODE=
QUIET_MODE=-q

HEURISTICS=
HEURISTICS=-D USE_BOUNDS
HEURISTICS=-D USE_REALCERTIFY

REALCERTIFY_PATH=~/Documents/GithubProjects/RealCertify
TEMP_LOG_FILE=${REALCERTIFY_PATH}/log_time.txt
SUFFIX=
OUTPUT_FILE=output${SUFFIX}.txt
LOG_FILE=log_time${SUFFIX}.txt
FILE=examples.mpl
FILE=benchmarking.mpl

QUIET_MODE=
HEURISTICS=-D USE_BOUNDS
FLAGS=-D LOG_TIME

.PHONY: all clean impl

all: impl

impl:
	if [ -f ${TEMP_LOG_FILE} ]; then rm ${TEMP_LOG_FILE}; fi;
	maple ${FILE} ${QUIET_MODE} ${HEURISTICS} ${FLAGS} > ${OUTPUT_FILE}
	if [ -f ${TEMP_LOG_FILE} ]; then mv ${TEMP_LOG_FILE} ${LOG_FILE}; fi;

clean:
	rm -rf ${OUTPUT_FILE}
	rm -rf ${LOG_FILE}
