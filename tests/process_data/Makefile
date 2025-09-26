PROG=./process_output.py
PROG=./process_time.py

PROG=./process_output.py

FILE=../data/experiments_8/output_archimedeanBenchmark.txt
FILE=../data/experiments_9/output_strictly_positive_polynomials_close_to_semialgebraic_sets.txt
FILE=../data/experiments_9/output_certificates_of_strictly_positive_polynomials.txt
FILE=../output_1_300.txt
FILE=../output_2_300.txt
FILE=../output_3_300.txt
FILE=../output_4_300.txt
FILE=../output_5_300.txt
FILE=../output.txt

.PHONY: certificates_of_strictly_positive_polynomials strictly_positive_polynomials_close_to_semialgebraic_sets

#all: certificates_of_strictly_positive_polynomials 
all: strictly_positive_polynomials_close_to_semialgebraic_sets

certificates_of_strictly_positive_polynomials:
	${PROG} ${FILE} 3 0

strictly_positive_polynomials_close_to_semialgebraic_sets:
	${PROG} ${FILE} 72 1
