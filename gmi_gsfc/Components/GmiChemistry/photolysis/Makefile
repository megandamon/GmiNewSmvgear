include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C include          all
	@$(MAKE) -C fastj            all
	@$(MAKE) -C fast_JX          all
	@$(MAKE) -C fast_JX65        all
	@$(MAKE) -C fast_JX53b       all
	@$(MAKE) -C fastJX53c_ref    all
	@$(MAKE) -C utils            all
	@$(MAKE) -C lookup           all
 
clean: 
	@$(MAKE) -C fastj            clean
	@$(MAKE) -C fast_JX          clean
	@$(MAKE) -C fast_JX65        clean
	@$(MAKE) -C fast_JX53b       clean
	@$(MAKE) -C fastJX53c_ref    clean
	@$(MAKE) -C utils            clean
	@$(MAKE) -C lookup           clean


#dependencies
