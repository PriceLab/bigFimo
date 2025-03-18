default:
	@echo roxy
	@echo install
	@echo build
	@echo test 
	@echo all [roxy install test]

roxy:
	- rm NAMESPACE
	R -e "devtools::document()"

install:
	R CMD INSTALL .  --no-test-load

build:
	R CMD build .
test:
	(cd inst/unitTests; make)

all: roxy install test


