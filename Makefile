PROJECT_NAME=StrictlyPositiveCert
SUFFIX=
FLAGS=

QUIET_MODE=
QUIET_MODE=-q

.PHONY: clean

all: ${PROJECT_NAME}.mla


${PROJECT_NAME}.mla: ${PROJECT_NAME}.mpl
	archive_maple_project.py ${PROJECT_NAME} ${PROJECT_NAME} ${PROJECT_NAME}

clean:
	rm -rf ${PROJECT_NAME}.mla
