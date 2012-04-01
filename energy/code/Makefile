OBJS = main.o function.o solver.o thomas.o
CC = g++
DEBUG = -g
LINK = -lm
CFLAGS = -Wall -c $(LINK) $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : energy main.o function.o solver.o thomas.o

energy : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o energy

main.o : main.cpp header.h
	$(CC) $(CFLAGS) main.cpp

function.o : function.cpp header.h
	$(CC) $(CFLAGS) function.cpp

solver.o : solver.cpp header.h
	$(CC) $(CFLAGS) solver.cpp

thomas.o : thomas.cpp header.h
	$(CC) $(CFLAGS) thomas.cpp

clean:
	\rm *.o energy

tar:
	tar cfvz energy2D.tgz main.cpp header.h function.cpp solver.cpp thomas.cpp
