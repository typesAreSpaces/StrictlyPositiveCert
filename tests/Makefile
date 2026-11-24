QUIET_MODE ?=
FLAGS ?=
SUFFIX ?=
FILE ?= benchmarking.mpl
OUTPUT_FILE ?= output$(SUFFIX).txt
LOG_FILE ?= log_time$(SUFFIX).txt

.PHONY: all library impl clean

all: impl

library:
	$(MAKE) -C ..

impl: library
	maple $(QUIET_MODE) $(FLAGS) "$(FILE)" > "$(OUTPUT_FILE)"

clean:
	rm -f -- "$(OUTPUT_FILE)" "$(LOG_FILE)"
