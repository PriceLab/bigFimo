default:
	@echo targets: roxy install test all

roxy:
	- rm NAMESPACE
	R -e "devtools::document()"

install:
	R CMD INSTALL .  --no-test-load
test:
	(cd inst/unitTests; make)

all: roxy install test


