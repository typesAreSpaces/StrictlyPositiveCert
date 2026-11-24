PYTHON ?= python3
PROG ?= ./process_isolated_output.py
FILE ?= ../output.txt

.PHONY: all certificates_of_strictly_positive_polynomials \
        strictly_positive_polynomials_close_to_semialgebraic_sets solving_stengles_example

all: certificates_of_strictly_positive_polynomials

certificates_of_strictly_positive_polynomials:
	$(PYTHON) "$(PROG)" "$(FILE)" 3 0

solving_stengles_example:
	$(PYTHON) "$(PROG)" "$(FILE)" 2 0

strictly_positive_polynomials_close_to_semialgebraic_sets:
	$(PYTHON) "$(PROG)" "$(FILE)" 72 1
