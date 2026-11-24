MAPLE ?= maple

FILE ?= to_discuss.mpl
QUIET_MODE ?= -q
# Use HEURISTIC=-DUSE_REALCERTIFY to select the alternate heuristic.
HEURISTIC ?= -D USE_BOUNDS
OUTPUT ?= output.txt

.PHONY: all clean

all:
	$(MAPLE) "$(FILE)" $(QUIET_MODE) $(HEURISTIC) > "$(OUTPUT)"

clean:
	rm -f -- "$(OUTPUT)"
