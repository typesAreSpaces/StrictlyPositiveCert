MAPLE ?= maple

# Examples: make FILE=examples.mpl QUIET_MODE=-q
# Use HEURISTICS=-DUSE_REALCERTIFY to select the alternate heuristic.
FILE ?= benchmarking.mpl
QUIET_MODE ?=
HEURISTICS ?= -D USE_BOUNDS
FLAGS ?= -D LOG_TIME
SUFFIX ?=
OUTPUT_FILE ?= output$(SUFFIX).txt
LOG_FILE ?= log_time$(SUFFIX).txt

.PHONY: all clean impl library

all: impl

library:
	$(MAKE) -C ..

impl: library
	$(MAPLE) "$(FILE)" $(QUIET_MODE) $(HEURISTICS) $(FLAGS) > "$(OUTPUT_FILE)"

clean:
	rm -f -- "$(OUTPUT_FILE)" "$(LOG_FILE)"
