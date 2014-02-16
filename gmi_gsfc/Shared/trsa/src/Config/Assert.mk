#                      -----------------                          
#                      Check Environment
#                      -----------------                          

ASSERT_FILE := $(wildcard $(ESMABIN)/Assert.pl)
ifeq ($(ASSERT_FILE),$(ESMABIN)/Assert.pl)
ifneq ($(shell $(ESMABIN)/Assert.pl), 0)
$(error Please correct your build environment and try again)
endif
else
$(warning WARNING: No environment check performed.)
endif
