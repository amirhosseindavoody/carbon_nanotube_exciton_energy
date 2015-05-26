
FC=gfortran
# FCFLAGS = -O3 -ffree-line-length-none
FCFLAGS += -fcheck=all -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination
FCFLAGS += -fbounds-check -ffree-line-length-none -fall-intrinsics
FCFLAGS += -I.

SRCDIR = ./src
OBJDIR = ./obj

A_exciton_energy_mod.o: comparams.o math_functions_mod.o write_log_mod.o
cnt_band_structure_mod.o: cnt_class.o comparams.o physical_constant_mod.o write_log_mod.o
cnt_geometry_mod.o: comparams.o math_functions_mod.o physical_constant_mod.o write_log_mod.o
comparams.o: cnt_class.o
dielectric_fn_mod.o: cnt_band_structure_mod.o comparams.o physical_constant_mod.o write_log_mod.o
E_exciton_energy_mod.o: comparams.o math_functions_mod.o write_log_mod.o
main.o: A_exciton_energy_mod.o cnt_band_structure_mod.o cnt_geometry_mod.o dielectric_fn_mod.o E_exciton_energy_mod.o parse_input_file_mod.o self_energy_mod.o write_log_mod.o
parse_input_file_mod.o: comparams.o physical_constant_mod.o write_log_mod.o
self_energy_mod.o: cnt_band_structure_mod.o comparams.o write_log_mod.o

main: main.o
	$(FC) $(FCFLAGS) -o $@.exe $(wildcard $(OBJDIR)/*.o) -llapack -lblas

%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) -c $(FCFLAGS) $< -J$(OBJDIR)
	@mv -f $@ $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.mod *.exe
	@rm -rf $(OBJDIR)