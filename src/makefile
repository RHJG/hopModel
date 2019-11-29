# specify the include directory (containing all .h files)
IDIR =../include
# specify compiler ... 
CC=gcc
# and compiler flags (-I$(IDIR): include the directory; -O3: pretty fast optimization).
CFLAGS=-I$(IDIR) -O3

# specify object directory to hide them from plain sight
ODIR=obj

# If the project depends on libraries specify them here
LIBS=-lm

# All headers 
_DEPS = functions.h main.h modelParameters.h variables.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# All objects that need to be created (for every .c-file we include a .o-file here)
_OBJ = functions.o initialize.o main.o memoryManager.o mersenneTwister.o outputManager.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Check for dependecies
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile the program
hoppingModel: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# If you need a fresh build call clean
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 