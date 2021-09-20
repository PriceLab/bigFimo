default:
	@echo targets: install test

install:
	R CMD INSTALL .
test:
	(cd inst/unitTests; make)


