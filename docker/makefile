IMG = bigfimobase
IMG2 = bigfimo2
DOCKER = sudo docker

default:
	@echo targets: build.1,2 bash.1,2 build2.noCache testClass test.rscript.in.container
build.1:
	$(DOCKER) build -t $(IMG) ./ -f Dockerfile.1

build.2:
	$(DOCKER) build -t $(IMG2) ./ -f Dockerfile.2

build.2.noCache:
	$(DOCKER) build --no-cache -t $(IMG2) ./ -f Dockerfile.2

bash.1:
	$(DOCKER) run -ti --rm $(IMG) bash
bash.2:
	$(DOCKER) run \
           -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results \
           -ti $(IMG2) bash \

testClass:
	sudo docker run -ti --rm bigfimo2 R -f /usr/local/scripts/test_BigFimoPackage.R

test.rscript.in.container:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 4040801 4041801 ZBTB7A 1 1e-5 \
	      /usr/local/data/human-jaspar2018-hocomoco-swissregulon.meme \
              /usr/local/data/results								

# a 2kb region around the tss, not yet in a container
outside.testScript:
	Rscript components/runBigFimo.R  chr1 161200670 161202670  NDUFS2 3 1e-6 \
             components/human-jaspar2018-hocomoco-swissregulon.meme \
        -v /data/dockerizedFlask/static/ampad/NIA-1898/liftover:/gatk/my_data              test.rscript.in.container


