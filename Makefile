# Compiler and flags
CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -g -pedantic-errors
LDFLAGS = -lm

C_FILES= symnmf.c

.PHONY: all clean

all: symnmf

symnmf:
	$(CC) $(CFLAGS) -o symnmf $(C_FILES) $(LDFLAGS)

clean:
	rm -f symnmf