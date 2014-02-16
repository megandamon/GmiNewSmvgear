include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include          all
	@$(MAKE) -C dao2utils        all
	@$(MAKE) -C dao2advec        all
	@$(MAKE) -C advectionMethod  all
	@$(MAKE) -C advecCoreMethod  all
 
clean: 
	@$(MAKE) -C dao2utils        clean
	@$(MAKE) -C dao2advec        clean
	@$(MAKE) -C advectionMethod  clean
	@$(MAKE) -C advecCoreMethod  clean


#dependencies
