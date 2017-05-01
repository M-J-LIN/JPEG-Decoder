all:main.c parser.c
	gcc -std=c99 -O2 parser.c main.c -o decoder -lm 
clean:
	rm decoder 
