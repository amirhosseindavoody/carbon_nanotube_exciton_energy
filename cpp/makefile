CC=g++-7
CC+= -O3

CFLAGS = -I./ -std=c++17
LFLAGS = -lm -lstdc++fs -std=c++17 -larmadillo

SRCDIR = ./src
OBJDIR = ./obj

object:
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $(SRCDIR)/*.cpp
	@mv -f ./*.o $(OBJDIR)

main: object
	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.exe
	@rm -rf $(OBJDIR)


#
#
# CC=g++
# CC+= -Ofast -finline-functions -ffp-contract=on
# # CC=icpc
# # CC+= -Ofast -simd -qopt-report-phase=vec
#
# CFLAGS = -I./ -I./lib/ -I./my_math/ -std=c++14
# # LFLAGS = -lgsl -lgslcblas -lm -std=c++14
# LFLAGS = -llapack -lblas -lm -std=c++14
#
# SRCDIR = ./src
# OBJDIR = ./obj
#
# object:
# 	@mkdir -p $(OBJDIR)
# 	$(CC) -c $(CFLAGS) $(SRCDIR)/*.cpp
# 	@mv -f ./*.o $(OBJDIR)
#
# main: object
# 	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)
#
# new: clean main
#
# # Utility targets
# .PHONY: clean
# clean:
# 	@reset
# 	@rm -f *.o *.exe
# 	@rm -rf $(OBJDIR)
