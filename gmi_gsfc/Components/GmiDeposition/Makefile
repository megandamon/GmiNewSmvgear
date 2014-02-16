include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include          all
	@$(MAKE) -C depositionMethod  all
 
clean: 
	@$(MAKE) -C depositionMethod  clean


#dependencies
