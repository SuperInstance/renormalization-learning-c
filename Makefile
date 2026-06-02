CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -Iinclude -O2 -g
LDFLAGS  = -lm

SRC      = src/lattice.c src/block_spin.c src/skill_model.c src/learning_curve.c
OBJ      = $(SRC:.c=.o)

TEST_SRC = tests/test_renorm.c
TEST_BIN = test_renorm

.PHONY: all test clean

all: $(OBJ)

$(OBJ): include/renorm_learn.h

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(OBJ) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(OBJ) $(TEST_BIN)
