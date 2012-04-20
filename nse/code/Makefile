OBJS = main.o cell.o field.o grid.o rect.o src.o thomas.o
CC = g++
DEBUG = -g
LINK = -lm
CFLAGS = -Wall -c $(LINK) $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : nse2d field.o cell.o rect.o grid.o src.o thomas.o main.o 

nse2d : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o nse2d

main.o : main.cc header.h constant.h
	$(CC) $(CFLAGS) main.cc

field.o : field.cc header.h constant.h
	$(CC) $(CFLAGS) field.cc

cell.o : cell.cc header.h constant.h
	$(CC) $(CFLAGS) cell.cc

rect.o : rect.cc header.h constant.h
	$(CC) $(CFLAGS) rect.cc

grid.o : grid.cc header.h constant.h
	$(CC) $(CFLAGS) grid.cc

src.o : src.cc header.h constant.h
	$(CC) $(CFLAGS) src.cc

thomas.o : thomas.c header.h constant.h
	$(CC) $(CFLAGS) thomas.c

clean:
	\rm *.o nse2d

tar:
	tar cfvz energy2D.tgz main.cc header.h field.cc cell.cc grid.cc src.cc thomas.c
