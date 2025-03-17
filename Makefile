CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -D_GNU_SOURCE

SRC_DIR = src
DOC_DIR = doc

OBJ = $(SRC_DIR)/main.o $(SRC_DIR)/functions.o

all: main run report.pdf

main: $(OBJ)
	$(CC) $(CFLAGS) -o main $(OBJ) -lm

$(SRC_DIR)/main.o: $(SRC_DIR)/main.c $(SRC_DIR)/functions.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/main.c -o $(SRC_DIR)/main.o

$(SRC_DIR)/functions.o: $(SRC_DIR)/functions.c $(SRC_DIR)/functions.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/functions.c -o $(SRC_DIR)/functions.o

run: main
	@echo "Uruchamiam program z parametrami: n_min=10, n_max=3000, step=100"
	@printf "10\n3000\n100\n" | ./main
	@echo "Program zakończył działanie."

report.pdf: $(DOC_DIR)/report.tex
	cd $(DOC_DIR) && pdflatex report.tex && pdflatex report.tex
	cp $(DOC_DIR)/report.pdf .

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(DOC_DIR)/*.aux $(DOC_DIR)/*.log $(DOC_DIR)/*.pdf $(DOC_DIR)/*.dat $(DOC_DIR)/*.out $(DOC_DIR)/*.png
	# rm -f main report.pdf
	
.PHONY: all run clean
