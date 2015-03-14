# this is simple make file

FC=gfortran
FCFLAGS = -O3 -ffree-line-length-none
#FCFLAGS = -Wall -Wtabs -Wconversion -fbounds-check -g -ffree-line-length-none
FCFLAGS += -I.

#this also works but I prefer the next set of commands
#main: module.f90 fnInput.f90 fnFileHandler.f90 smallFns.f90 fnPhysConst.f90 fnGeomProp.f90 fnCNTband.f90 fnGrapheneEnergy.f90 fnDielectric.f90 fnSelfEnergy.f90 fnExcitonDispersion.f90 fnExcitonEnergy.f90 main.f90

main.o: main.f90 module.o fnInput.o smallFns.o fnPhysConst.o fnGeomProp.o fnCNTband.o fnGrapheneEnergy.o fnDielectric.o fnSelfEnergy.o fnExcitonDispersion.o fnExcitonEnergy.o
main: module.o fnInput.o fnFileHandler.o smallFns.o fnPhysConst.o fnGeomProp.o fnCNTband.o fnGrapheneEnergy.o fnDielectric.o fnSelfEnergy.o fnExcitonDispersion.o fnExcitonEnergy.o
fnInput.o: fnInput.f90
fnFileHandler.o : fnFileHandler.f90
smallFns.o : smallFns.f90
fnPhysConst.o : fnPhysConst.f90
fnGeomProp.o : fnGeomProp.f90
fnCNTband.o : fnCNTband.f90
fnGrapheneEnergy.o : fnGrapheneEnergy.f90
fnDielectric.o : fnDielectric.f90
fnSelfEnergy.o : fnSelfEnergy.f90
fnExcitonDispersion.o : fnExcitonDispersion.f90
fnExcitonEnergy.o : fnExcitonEnergy.f90

%: %.f90
	$(FC) -o $@ $^ $(FCFLAGS) -llapack -lblas

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *~ main