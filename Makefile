PROJECT_NAME=StrictlyPositiveCert

# OPTIONS:
#FLAGS+=-DENABLE_VERIFICATION
#FLAGS+=-DENABLE_BINARY_SEARCH_AVKL
#FLAGS+=-DENABLE_AVERKOV_CHECK
#FLAGS+=-DENABLE_QUICK_CERT
#FLAGS+=-DLOG_TIME

#DEBUG_MODE+=-DENABLE_DEBUGGING

#VALUES+=-DEPS_FACTOR=<rational>
#VALUES+=-DN_GUESS_AVKL=<integer>

FLAGS=
DEBUG_MODE=
VALUES=

FLAGS+=-DENABLE_BINARY_SEARCH_AVKL
VALUES+=-DEPS_FACTOR=17/20
VALUES+=-DN_GUESS_AVKL=70
DEBUG_MODE+=-DENABLE_DEBUGGING

.PHONY: clean test lint

all: ${PROJECT_NAME}.mla

${PROJECT_NAME}.mla: $(wildcard src/*.mpl) scripts/installation_script.mpl
	maple -q ${FLAGS} ${DEBUG_MODE} ${VALUES} ./scripts/installation_script.mpl

clean:
	rm -rf ${PROJECT_NAME}.mla

lint:
	find src/ -type f -name '*.mpl' -print0 | xargs -0 -n 1 format-maple.sh

# Maple can exit successfully after reporting an error; the runner checks output too.
test:
	python3 tests/run_regressions.py ${FLAGS} ${VALUES}
