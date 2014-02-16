include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include          all
	@$(MAKE) -C sad              all
	@$(MAKE) -C mechanisms       all
	@$(MAKE) -C solvers          all
	@$(MAKE) -C photolysis       all
	@$(MAKE) -C AerosolDust      all
	@$(MAKE) -C ioChemistry      all
	@$(MAKE) -C chemistryMethod  all
 
clean: 
	@$(MAKE) -C sad              clean
	@$(MAKE) -C mechanisms       clean
	@$(MAKE) -C solvers          clean
	@$(MAKE) -C photolysis       clean
	@$(MAKE) -C AerosolDust      clean
	@$(MAKE) -C ioChemistry      clean
	@$(MAKE) -C chemistryMethod  clean


#dependencies
