CC=g++
CC+= -Ofast -finline-functions -ffp-contract=on
# CC=icpc
# CC+= -Ofast -simd -qopt-report-phase=vec

CFLAGS = -I./ -I./lib/ -I./my_math/ -std=c++14
# LFLAGS = -lgsl -lgslcblas -lm -std=c++14
LFLAGS = -llapack -lblas -lm -std=c++14

SRCDIR = ./src
OBJDIR = ./obj

object:
	@mkdir -p $(OBJDIR)
	$(CC) -c $(CFLAGS) $(SRCDIR)/*.cpp
	@mv -f ./*.o $(OBJDIR)

main: object
	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)

new: clean main

# Utility targets
.PHONY: clean
clean:
	@reset
	@rm -f *.o *.exe
	@rm -rf $(OBJDIR)