CC = gcc
CFLAGS = -Igraphics -Wall -Wextra -O3 #-march=native -ffast-math -funroll-loops -flto -ftree-vectorize -g
TARGET = galsim
LDLIBS = -lX11 -lm -fopenmp 

all: $(SUB_TARGET) $(TARGET)

$(TARGET): galsim.c $(SUB_TARGET)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) galsim.o $(SUB_TARGET)
	@for dir in $(SUB_DIR); do $(MAKE) -C $$dir clean; done

.PHONY: all clean
