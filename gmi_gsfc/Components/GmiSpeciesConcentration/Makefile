include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C ioSpcConcentration      all
	@$(MAKE) -C spcConcentrationMethod  all
 
clean: 
	@$(MAKE) -C ioSpcConcentration      clean
	@$(MAKE) -C spcConcentrationMethod  clean


#dependencies
