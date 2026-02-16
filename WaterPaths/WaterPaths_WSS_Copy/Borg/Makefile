# CC = gcc
# MPICC = mpicc
# CFLAGS = -O3
# LDFLAGS = -Wl,-R,\.
# LIBS = -lm
# UNAME_S = $(shell uname -s)

# ifneq (, $(findstring SunOS, $(UNAME_S)))
# 	LIBS += -lnsl -lsocket -lresolv
# else ifneq (, $(findstring MINGW, $(UNAME_S)))
# 	# MinGW is not POSIX compliant
# else
# 	POSIX = yes
# endif

# compile:
# 	$(CC) $(CFLAGS) $(LDFLAGS) -o dtlz2_serial.exe dtlz2_serial.c borg.c mt19937ar.c $(LIBS)

# ifdef POSIX
# 	$(CC) $(CFLAGS) $(LDFLAGS) -o borg.exe frontend.c borg.c mt19937ar.c $(LIBS)
# endif

# 	$(MPICC) $(CFLAGS) $(LDFLAGS) -o dtlz2_ms.exe dtlz2_ms.c borgms.c mt19937ar.c $(LIBS)
# 	$(MPICC) $(CFLAGS) $(LDFLAGS) -o dtlz2_mm.exe dtlz2_mm.c borgmm.c mt19937ar.c $(LIBS)
# .PHONY: compile

CC = gcc
ARGS = -Wall -O3 
MPICC = mpicc
MPIARGS = -O3
CFLAGS = -O3
LDFLAGS = -Wl,-R,\.

all:
	${CC} ${ARGS} -ansi -o mt19937ar.o -c mt19937ar.c
	${CC} ${ARGS} -ansi -o borg.o -c borg.c
	ar rcs libborg.a borg.o mt19937ar.o
	${MPICC} ${MPIARGS} -o borgmm.o -c borgmm.c
	ar rcs libborgmm.a borgmm.o mt19937ar.o
	${MPICC} ${MPIARGS} -o borgms.o -c borgms.c
	ar rcs libborgms.a borgms.o mt19937ar.o

clean: 
	rm *.o *.a *.exe *.so

