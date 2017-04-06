CC=g++
CC+= -Wno-deprecated-declarations

CFLAGS = -I./ -I../ -std=c++14
LFLAGS = -llapack -lblas -lm -std=c++14

SRCDIR = .
OBJDIR = .

object:
	$(CC) $(CFLAGS) -c $(SRCDIR)/*.cpp

test_eigen: object
	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)

# Utility targets
.PHONY: clean
clean:
	@reset
	@rm -f *.o *.exe
	# @rm -rf $(OBJDIR)