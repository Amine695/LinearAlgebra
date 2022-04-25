CC=gcc
CFLAGS= -g -Wall -pedantic -O3
LDFLAGS = -lm
SRC= $(wildcard src/*.c)  
BUILD = $(build)
OBJ= $(SRC:.c=.o)
EXEC = project

all: $(EXEC)
(build)/%.o: %.c   
	$(CC) $(CFLAGS) -o $@ -c $<   
	 
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@   
	

clean: 
	rm -f $(OBJ) $(EXEC) 	

