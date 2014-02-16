include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all: 
	@$(MAKE) -C sharedVariables  all
	@$(MAKE) -C $(SULFCASE)      all
	@$(MAKE) -C smv2chem         all
	@$(MAKE) -C ammonia          all
 
clean: 
	@$(MAKE) -C sharedVariables  clean
	@$(MAKE) -C $(SULFCASE)      clean
	@$(MAKE) -C smv2chem         clean
	@$(MAKE) -C ammonia          clean


#dependencies
