ifdef $(ARCHi) # user has defined the architecture
else           # or determine the architecture automatically
  ARCHi = $(shell uname)

  ChemMecha = $CHEMCASE
endif

machineName = $(shell hostname | awk '{print substr($0,0,8)}')

#export machineName
