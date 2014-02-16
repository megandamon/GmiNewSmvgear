include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C GmiSpeciesConcentration all
	@$(MAKE) -C GmiEmission             all
	@$(MAKE) -C GmiDeposition           all
	@$(MAKE) -C GmiConvection           all
	@$(MAKE) -C GmiDiffusion            all
	@$(MAKE) -C GmiAdvection            all
	@$(MAKE) -C GmiChemistry            all
 

clean: 
	@$(MAKE) -C GmiSpeciesConcentration clean
	@$(MAKE) -C GmiEmission             clean
	@$(MAKE) -C GmiDeposition           clean
	@$(MAKE) -C GmiConvection           clean
	@$(MAKE) -C GmiDiffusion            clean
	@$(MAKE) -C GmiAdvection            clean
	@$(MAKE) -C GmiChemistry            clean


#dependencies
