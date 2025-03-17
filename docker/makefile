IMG = bigfimobase
IMG2 = pshannon/bigfimo
DOCKER = sudo docker

default:
	@echo build.1
	@echo build.2
	@echo bash.1
	@echo bash.2 
	@echo build2.noCache 
	@echo testClass
	@echo test.rscript.in.container.small
	@echo test.rscript.in.container.medium
	@echo test.rscript.in.container.huge
	@echo run.klf1.2k.noDocker
	@echo run.klf1.2k
	@echo run.klf1.20k
	@echo run.klf1.100k
	@echo run.klf1.500k
	@echo run.klf1.1000k

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
           -v /users/pshannon/github/bigFimo/docker/inputs:/usr/local/data/inputs \
           -ti $(IMG2) bash


push:
	$(DOCKER) push $(IMG2)

testClass:
	sudo docker run -ti --rm bigfimo2 R -f /usr/local/scripts/test_BigFimoPackage.R

# 8:02 1000 bases, 1e-5, less than 60 seconds
test.script.in.container.small:
	$(DOCKER) run \
          -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results \
          -v /users/pshannon/github/bigFimo/docker/inputs:/usr/local/data/inputs \
          -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 4040801 4041801 ZBTB7A 2 1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results	

# 2 processes, over 1kb, finds 62 hits, ~2kb file
test.rscript.in.container.small:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 4040801 4041801 ZBTB7A 2 1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

# 5 processes, 68kb: 715 fimo hits saved in a file ~20k 
test.rscript.in.container.medium:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 161165894 161234221 NDUFS2 5 1e-6 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

# 20 processes, 1093kb:  2,737,622  fimo hits save in 46,443,187 file, 12 minutes
test.rscript.in.container.huge:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 160653434 161746681 NDUFS2 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

test.rscript.withChromLocTable:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimoWithChromLocsTable.R \
              /usr/local/data/pilraMicro.csv \
              PILRA \
              5 \
              1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

# tss is chr19 12,887,201
# 12888201 - 12886201 = 2000 bases only
run.klf1.2k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19  12886201 12888201 KLF1 3 1e-5 \
	      /usr/local/data/human-jaspar2024-hocomocov13.meme \
              /usr/local/data/results								
# 	/usr/local/data/human-jaspar2022-hocomoco-core-A.meme \


run.klf1.2k.noDocker:
	Rscript components/runBigFimo.R chr19  12886201 12888201 \
          KLF1 3 1e-5 \
          components/human-jaspar2024-hocomocov13.meme results/

# ~1 minute, 118086 rows
run.klf1.20k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 12867201 12907201 KLF1 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results

# 2 minutes, 292144 rows
# tss chr19 12,887,201
# 
run.klf1.100k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 12837201 12937201  KLF1 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results

# minutes: 7:28
# rows: 1576717
# tss chr19 12,887,201
# 
run.klf1.500k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 12637201 13137201  KLF1 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results

# minutes: 
# rows: 
# tss chr19 12,887,201
# 
run.klf1.1000k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 1576717 1576717 KLF1 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results

# chr1:47,216,290-47,232,335
run.tal1.2k:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 47232335 47234335 TAL1 20 1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results


# chr1:47,216,290-47,232,335
run.tal1.2M:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 46232335 48232335 TAL1 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results

