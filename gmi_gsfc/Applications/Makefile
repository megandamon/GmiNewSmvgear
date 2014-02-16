include $(CONFIG)/libraries.mk
include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all:
	@$(MAKE) -C GmiOutputs  all
	@$(MAKE) -C GmiUtils    all
	@$(MAKE) -C GmiMain     gmi.x


clean:
	@$(MAKE) -C GmiOutputs  clean
	@$(MAKE) -C GmiUtils    clean
	@$(MAKE) -C GmiMain     clean
	rm -f GmiBin/*.x
