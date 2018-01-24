CC=g++-7
CC+= -O3
# CC+= -Wall -Wno-comment -ansi -pedantic-errors -g

CFLAGS = -I./ -std=c++17
LFLAGS = -lstdc++fs -std=c++17 -larmadillo

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