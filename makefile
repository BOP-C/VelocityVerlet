CC=g++
CFLAGS=-march=native -O2
OBJECTS = velocityVerlet.o point.o position.o velocity.o
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
TARGET = velocityVerlet
else
TARGET = velocityVerlet.exe
endif

velocityVerlet : $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

$(OBJECTS) : common.h
velocityVerlet.o point.o : point.h
position.o : position.h
velocity.o : velocity.h

.PHONY : clean
clean :
	-rm $(TARGET) $(OBJECTS)
