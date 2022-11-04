IMG = bigfimo:03nov2022
IMG2 = bigfimo2:03nov2022
DOCKER = sudo docker

default:
	@echo targets: build.[1|2] bash.[1|2]

build.1:
	$(DOCKER) build -t $(IMG) ./ -f Dockerfile.1

build.2:
	$(DOCKER) build -t $(IMG2) ./ -f Dockerfile.2

buildNoCache:
	$(DOCKER) build --no-cache -t $(IMG) ./

bash.1:
	$(DOCKER) run -ti --rm $(IMG) bash
bash.2:
	$(DOCKER) run -ti --rm $(IMG2) bash

host.install:
	R -f installTrenaPrequisites.R

rstudio.bioc:
	$(DOCKER) run -e PASSWORD=password -p 8787:8787 bioconductor/release_core2

rstudio.trena:
	$(DOCKER) run -e PASSWORD=password -p 8787:8787 $(ING)

