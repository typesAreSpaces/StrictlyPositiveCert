PROJECT_NAME ?= StrictlyPositiveCert
MAPLE ?= maple
PYTHON ?= python3

# Override these on the command line to customize the library build.
# Optional flags: -DENABLE_VERIFICATION, -DENABLE_AVERKOV_CHECK,
# -DENABLE_QUICK_CERT, -DLOG_TIME.
FLAGS ?= -DENABLE_BINARY_SEARCH_AVKL
DEBUG_MODE ?= -DENABLE_DEBUGGING
VALUES ?= -DN_GUESS_AVKL=70

.PHONY: all clean lint test

all: $(PROJECT_NAME).mla

$(PROJECT_NAME).mla: $(wildcard src/*.mpl) scripts/installation_script.mpl
	$(MAPLE) -q $(FLAGS) $(DEBUG_MODE) $(VALUES) ./scripts/installation_script.mpl

clean:
	rm -f -- "$(PROJECT_NAME).mla"

lint:
	find src/ -type f -name '*.mpl' -print0 | xargs -0 -n 1 format-maple.sh

# Maple can exit successfully after reporting an error; the runner checks output too.
test:
	$(PYTHON) tests/run_regressions.py $(FLAGS) $(VALUES)
