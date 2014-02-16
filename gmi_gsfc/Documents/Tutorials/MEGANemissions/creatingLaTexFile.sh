#
source ~/.cshrc.gmi
#
rm -f gmiMEGANemissions*
#
./protex -s anIntroduction $gmi/Components/GmiEmission/MEGAN/GmiEmissionMEGAN_mod.F90 $gmi/Components/GmiEmission/ioEmission/ReadInputMEGAN_mod.F90 $gmi/Shared/GmiMetFields/GmiSurfaceTemperature_mod.F90 > gmiMEGANemissions.tex
#

