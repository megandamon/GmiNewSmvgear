include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include_setkin   all
	@$(MAKE) -C setkin           all

clean: 
	@$(MAKE) -C setkin           clean


#dependencies
