CC=g++
INCLUDE=-I/usr/include/GL -I/home/mos/linux_Programming/lib/ann_1.1.2/ann_1.1.2/include
LIB=-L/usr/X11R6/lib -L/usr/lib
OUTBIN=line

all: clean main.o
#	${CC} ${LIB} *.o -lGL -lglut -lGLU -lstdc++ -lm -lann -o ${OUTBIN}
	${CC} ${LIB} *.o -lGL -lglut -lGLU -lstdc++ -lm -o ${OUTBIN}

main.o:
	${CC} ${INCLUDE} *.cpp -c -Wunused -Wno-deprecated -g

clean:
	rm -rf *.o
	rm -rf *.*~
	rm -rf ${OUTBIN}
