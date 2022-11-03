IMG = bigfimo:03nov2022

default:
	@echo targets: build bash 

build:
	docker build -t $(IMG) ./

buildNoCache:
	docker build --no-cache -t $(IMG) ./

bash:
	docker run -ti --rm $(IMG) bash

host.install:
	R -f installTrenaPrequisites.R

rstudio.bioc:
	docker run -e PASSWORD=password -p 8787:8787 bioconductor/release_core2

rstudio.trena:
	docker run -e PASSWORD=password -p 8787:8787 $(ING)

