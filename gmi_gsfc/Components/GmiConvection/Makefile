include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C convectionMethod  all
 
clean: 
	@$(MAKE) -C convectionMethod  clean


#dependencies
