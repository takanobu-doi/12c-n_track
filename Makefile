CXX = g++
TARGET = my_track
SOURCE = main.o mktrack.o gen_eve.o database.o dataset.o nuclear.o
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
GARFIELDFLAGS = -I${GARFIELD_HOME}/Include
GARFIELDLIBS = ${GARFIELD_HOME}/Library/libGarfield.a

CFLAGS = -O2 -Wall ${ROOTFLAGS} ${GARFIELDFLAGS}
LIBS = ${GARFIELDLIBS} -lgfortran ${ROOTLIBS}
DEBAG = -g

all: ${TARGET}
${TARGET}: ${SOURCE}
	${CXX} $^ -o $@ ${CFLAGS} ${LIBS}
.cpp.o:
	${CXX} -c ${CFLAGS} $<
clean:
	${RM} *.o ${TARGET} *~
reflesh:
	${RM} *.o *~
