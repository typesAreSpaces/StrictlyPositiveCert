QUIET_MODE ?= -q
FLAGS ?=
FILE ?= to_discuss.mpl
OUTPUT ?= output.txt

.PHONY: all clean

all:
	maple $(QUIET_MODE) $(FLAGS) "$(FILE)" > "$(OUTPUT)"

clean:
	rm -f -- "$(OUTPUT)"
