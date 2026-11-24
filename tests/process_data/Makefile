# Override PROG or FILE to process a different benchmark output.
PROG ?= ./process_isolated_output.py
FILE ?= ../output.txt

.PHONY: all certificates_of_strictly_positive_polynomials solving_stengles_example strictly_positive_polynomials_close_to_semialgebraic_sets

all: certificates_of_strictly_positive_polynomials

certificates_of_strictly_positive_polynomials:
	$(PROG) "$(FILE)" 3 0

solving_stengles_example:
	$(PROG) "$(FILE)" 2 0

strictly_positive_polynomials_close_to_semialgebraic_sets:
	$(PROG) "$(FILE)" 72 1
