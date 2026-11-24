PROJECT_NAME := StrictlyPositiveCert

# Optional FLAGS: -DENABLE_VERIFICATION, -DENABLE_AVERKOV_CHECK,
# -DENABLE_QUICK_CERT, -DLOG_TIME. Set FLAGS= to disable binary search.
FLAGS ?= -DENABLE_BINARY_SEARCH_AVKL
DEBUG_MODE ?= -DENABLE_DEBUGGING
VALUES ?= -DN_GUESS_AVKL=70

SOURCES := $(wildcard src/*.mpl)
INSTALL_SCRIPT := scripts/installation_script.mpl

.PHONY: all clean test lint

all: $(PROJECT_NAME).mla

$(PROJECT_NAME).mla: $(SOURCES) $(INSTALL_SCRIPT)
	maple -q $(FLAGS) $(DEBUG_MODE) $(VALUES) "$(INSTALL_SCRIPT)"

clean:
	rm -f -- "$(PROJECT_NAME).mla"

lint:
	find src/ -type f -name '*.mpl' -print0 | xargs -0 -n 1 format-maple.sh

# Maple can exit successfully after reporting an error; the runner checks output too.
test:
	python3 tests/run_regressions.py $(FLAGS) $(VALUES)
