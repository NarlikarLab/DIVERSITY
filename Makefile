CC = gcc -g -O3
# CFLAGS1 = -w
CFLAGS1 = -Wall -Wextra -ansi -pedantic-errors
CFLAGS = -fPIC -Wall 
LDFLAGS = -shared
LM = -lm
RM = rm -f
TARGET_LIB = libdatacal.so
QSRCS = messages.c linkedListOps.c modelops.c motifops.c traindata.c
VAR = -L$(CURDIR) -Wl,-rpath=$(CURDIR)
# VAR = -rpath=$(CURDIR)
OBJS = $(QSRCS:.c=.o)

.PHONY: all

all: ${TARGET_LIB}
	cat cst1 > cstructures.py; echo "libctest = cdll.LoadLibrary(\"$(CURDIR)/libdatacal.so\")" >> cstructures.py; cat cst2 >> cstructures.py
	echo "python $(CURDIR)/diversityMain.py "\"\$$@\"" $(CURDIR)/" > diversity; chmod +x diversity

$(TARGET_LIB): $(OBJS) 
	   $(CC) ${LDFLAGS} -o $@ $^ -lm

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(QSRCS:.c=.d) $(QFILE) *~ *.pyc diversity cstructures.py 
