CC	= gcc
BINDIR	= ../bin
INCDIR	= ../inc

CFLAGS	= -O -I${INCDIR}
CLIBS	= -lnsl -lm
#CLIBS	= -lm
OBJS	= main.o HF_EWO_fft.o

TARGET	= ${BINDIR}/main

${TARGET}: ${OBJS}
	${CC} ${CFLAGS} -o $@ ${OBJS} ${CLIBS}

.c.o:
	${CC} ${CFLAGS} -c $< ${CLIBS} 

clean:
	rm -f *.o *~ ${TARGET}
