include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C $(SETKINDIR)     all
 
clean: 
	@$(MAKE) -C $(SETKINDIR)     clean


#dependencies
