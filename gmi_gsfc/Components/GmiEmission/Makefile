include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include           all
	@$(MAKE) -C diagnEmission     all
	@$(MAKE) -C GOCARTroutines    all
	@$(MAKE) -C ioEmission        all
	@$(MAKE) -C MEGAN             all
	@$(MAKE) -C GalacticCosmicRay all
	@$(MAKE) -C llnl              all
	@$(MAKE) -C harvard           all
	@$(MAKE) -C lightning         all
	@$(MAKE) -C emissionMethod    all
 
clean: 
	@$(MAKE) -C GOCARTroutines    clean
	@$(MAKE) -C GalacticCosmicRay clean
	@$(MAKE) -C llnl              clean
	@$(MAKE) -C harvard           clean
	@$(MAKE) -C lightning         clean
	@$(MAKE) -C ioEmission        clean
	@$(MAKE) -C diagnEmission     clean
	@$(MAKE) -C MEGAN             clean
	@$(MAKE) -C emissionMethod    clean

#dependencies
