include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C GmiInclude           all
	@$(MAKE) -C GmiIOutilities       all
	@$(MAKE) -C GmiESMF              all
	@$(MAKE) -C NcUtils_Single       all
	@$(MAKE) -C GmiCommunications    all
	@$(MAKE) -C GmiDomainDecomp      all
	@$(MAKE) -C GmiSupportingModules all
	@$(MAKE) -C GmiDiagnostics       all
	@$(MAKE) -C GmiMetFields         all
 

clean: 
	@$(MAKE) -C GmiIOutilities       clean
	@$(MAKE) -C GmiESMF              clean
	@$(MAKE) -C NcUtils_Single       clean
	@$(MAKE) -C GmiCommunications    clean
	@$(MAKE) -C GmiDomainDecomp      clean
	@$(MAKE) -C GmiSupportingModules clean
	@$(MAKE) -C GmiDiagnostics       clean
	@$(MAKE) -C GmiMetFields         clean


#dependencies
