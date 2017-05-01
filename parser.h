#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define SOI_SYMB 0xffd8
#define EOI_SYMB 0xffd9
#define APP_SYMB 0xffe0
#define DQT_SYMB 0xffdb
#define SOF_SYMB 0xffc0
#define DHT_SYMB 0xffc4
#define SOS_SYMB 0xffda
const static float cu[8] = {0.707107, 1, 1, 1, 1, 1, 1, 1};
const static float cv[8] = {0.707107, 1, 1, 1, 1, 1, 1, 1};
const static float cos_tb[8][8] = {
  1.000000,   0.980785,   0.923880,   0.831470,   0.707107,   0.555570,   0.382683,   0.195090, 
  1.000000,   0.831470,   0.382683,  -0.195090,  -0.707107,  -0.980785,  -0.923880,  -0.555570,
  1.000000,   0.555570,  -0.382683,  -0.980785,  -0.707107,   0.195090,   0.923880,   0.831470,
  1.000000,   0.195090,  -0.923880,  -0.555570,   0.707107,   0.831470,  -0.382683,  -0.980785, 
  1.000000,  -0.195090,  -0.923880,   0.555570,   0.707107,  -0.831470,  -0.382683,   0.980785, 
  1.000000,  -0.555570,  -0.382683,   0.980785,  -0.707107,  -0.195090,   0.923880,  -0.831470, 
  1.000000,  -0.831470,   0.382683,   0.195090,  -0.707107,   0.980785,  -0.923880,   0.555570, 
  1.000000,  -0.980785,   0.923880,  -0.831470,   0.707107,  -0.555570,   0.382683,  -0.195090 
};
typedef struct APP{
	int length;
}APP;

typedef struct DQT{
	int length;
	int pq, tq; // pq:precision(0->Qk=8 ~ 1->Qk=16), tq:ID(0~3)
	int* table;
	int* zz;
}DQT;
typedef struct NF{
	int c; // component size
	int h, v, tq; // h:horizontal sampling , v:vertical sampling , tq:quantization table ID
}NF;
typedef struct SOF{
	int length;
	int p; // precision
	int y, x; // y:high, x:width
	int nfSize;
	NF *nf;
}SOF;

typedef struct DHT{
	int length;
	int tc, th; // tc:0->DC;1->AC , th:huffman table ID
	uint8_t li[16];
	uint8_t accu[16];
	uint8_t* vi; // content
	int size;
	uint16_t* codeword;
	uint16_t* codelength;
}DHT;

typedef struct NS{
	uint8_t cs, td, ta;
}NS;

typedef struct SOS{
	int length;
	int nsSize; 
	NS* ns;
	uint8_t ss, se, ah, al;
}SOS;
typedef struct YCbCr{
	int table[4][64];
	int zz[4][8][8];
	int idct[4][8][8];
}YCbCr;

void parse(char * filename);
