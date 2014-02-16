#
source ~/.cshrc.gmi
#
cp $gmi/Shared/GmiCommunications/*F90 .
#
./protex -s anIntroduction GmiDomainDecomposition_mod.F90 GmiGrid_mod.F90 GmiSubdomain2Procs_mod.F90 GmiSub2Glob_mod.F90 GmiGlob2Sub_mod.F90 GmiReduceSlaves_mod.F90 GmiBroadcast_mod.F90 GmiSubDomainsBC_mod.F90 GmiGhostZones_mod.F90 GmiReduce_mod.F90 GmiMessageInit_mod.F90 GmiMessageFinal_mod.F90 GmiRootProcessor_mod.F90 GmiMessagePassing_mod.F90 > removalESMpackage.tex
#
rm -f *.F90

