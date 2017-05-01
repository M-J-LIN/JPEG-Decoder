#include <stdio.h>
#include <stdlib.h>
#include "parser.h"

int main(int argc, char* argv[]){
	if(argc != 2){
		printf("need an input jpg file\n");             
  		return 0;
    }
	parse(argv[1]);

	return 0;
}

