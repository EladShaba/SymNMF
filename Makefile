# Compiler and flags
C = gcc
FLAG = -ansi -Wall -Wextra -Werror -g -pedantic-errors
LF = -lm

C_FILES= symnmf.c

.PHONY: all clean

all: symnmf

symnmf:
	$(C) $(FLAG) -o symnmf $(C_FILES) $(LF)

