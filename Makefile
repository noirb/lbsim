# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -pedantic -Wall -std=c99
DBGFLAGS=-g3 -gdwarf-3

# Linker flags
# ------------
LDFLAGS= -lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: CFLAGS += -DNDEBUG
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

debug: CFLAGS += $(DBGFLAGS)
debug: $(EXECUTABLE)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
