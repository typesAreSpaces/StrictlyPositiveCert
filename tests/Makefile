QUIET_MODE=
QUIET_MODE=-q

HEURISTICS=
HEURISTICS=-D USE_BOUNDS
HEURISTICS=-D USE_REALCERTIFY

SUFFIX=
OUTPUT_FILE=output${SUFFIX}.txt
LOG_FILE=log_time${SUFFIX}.txt
FILE=benchmarking.mpl
FILE=examples.mpl
FILE=comparison_benchmark.mpl

QUIET_MODE=
HEURISTICS=-D USE_BOUNDS
FLAGS=-D LOG_TIME

.PHONY: all clean impl

all: impl

impl:
	maple ${FILE} ${QUIET_MODE} ${HEURISTICS} ${FLAGS} > ${OUTPUT_FILE}

clean:
	rm -rf ${OUTPUT_FILE}
	rm -rf ${LOG_FILE}
