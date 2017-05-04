#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "parser.h"

static FILE* fptr;
static uint8_t bmp[1024][2048][3];
static uint8_t buffer[2];
static char outputS[50];
APP app;
DQT dqt[4];
int dqt_count = 0;
SOF sof;
DHT dht[4];
int dht_count = 0;
SOS sos;
int maxY, maxX;
int max(int a, int b){return (a>b) ? a : b;}
void parseAPP(){
	//printf("%x%x - APP\n", buffer[0], buffer[1]);
	if(fread(buffer, 1, 2, fptr)){
		app.length = buffer[0]*256 + buffer[1];
		//printf("        Lq: %d\n", app.length);
	}
	fseek ( fptr , app.length - 2 , SEEK_CUR );
}
void zigzag(int* dest, int* source, int max){
	dest[0] = source[0];
	int a = 0, b = 1, da = 1, db = -1, count = 1;
	for(int i = 1; i < 14; i++){
		int num;
		if(i >= 8)
			num = 14-i, a += da, b +=db;
		else 
			num = i;
		for(int j = 0; j <= num; j++){
			if(count >= max)
				return;
			dest[a * 8 + b] = source[count];
			a += da, b += db, count++;
		}
		if(da == 1)
			b++;
		else
			a++;
		da *= -1, db *= -1;
	}
	dest[63] = source[63];
	
}
void parseDQT(int id){
	//printf("%x%x - DQT\n", buffer[0], buffer[1]);
	int length;
	if(fread(buffer, 1, 2, fptr)){
		dqt[id].length = buffer[0]*256 + buffer[1];
		length = dqt[id].length;
		//printf("        Lq: %d\n", dqt[id].length);
	}
	for(int i = id; i < id + (length - 2) / 65; i++){
		if(fread(buffer, 1, 1, fptr)){
			dqt[i].pq = buffer[0] >> 4, dqt[i].tq = buffer[0] % 16;
			//printf("        Pq: %d    Tq: %d    Qk:8\n", dqt[i].pq, dqt[i].tq);
		}
		//printf("        In Parsing order (Zig-Zag-Scan order):\n          ");
		dqt[i].table = (int *)malloc(sizeof(int) * dqt[0].length);
		dqt[i].zz = (int *)malloc(sizeof(int) * dqt[0].length);
		memset((void *)dqt[i].zz, 0, sizeof(dqt[i].zz));
		for(int j = 0; j < 64; j++){
			uint8_t tmp;
			fread(&tmp, 1, 1, fptr); 
			dqt[i].table[j] = tmp;
			//printf("%d ",dqt[i].table[j]);
		}
		//printf("\n");
		//printf("        In Block order:\n          ");
		zigzag(dqt[i].zz, dqt[i].table, 64);
		for(int j = 0; j < 8; j++){
			for(int k = 0; k < 8; k++){
				//printf("%2d ", dqt[i].zz[j * 8 + k]);
			}
			//printf("\n          ");
		}
		//printf("\n");

	}
}
void parseSOF(){
	//printf("%x%x - SOF\n", buffer[0], buffer[1]);
	if(fread(buffer, 1, 2, fptr)){
		sof.length = buffer[0]*256 + buffer[1];
		//printf("        Lf: %d", sof.length);
	}
	if(fread(buffer, 1, 1, fptr)){
		sof.p = buffer[0];
		//printf("  P: %d", sof.p);
	}
	if(fread(buffer, 2, 1, fptr)){
		sof.y = buffer[0]*256+buffer[1];
		//printf("  Y: %d", sof.y);
	}
	if(fread(buffer, 2, 1, fptr)){
		sof.x = buffer[0]*256+buffer[1];
		//printf("  X: %d", sof.x);
	}
	if(fread(buffer, 1, 1, fptr)){
		sof.nfSize = buffer[0];
		//printf("  Nf: %d\n", sof.nfSize);
	}
	sof.nf = (NF *)malloc(sizeof(NF) * sof.nfSize);
	uint8_t tmp[3 * sof.nfSize];
	fread(tmp, 3 * sof.nfSize, 1, fptr);
	for(int i = 0; i < sof.nfSize; i++){
		sof.nf[i].c = tmp[i * 3], sof.nf[i].h = tmp[i * 3 + 1] >> 4, 
		sof.nf[i].v = tmp[i * 3 + 1] % 16, sof.nf[i].tq = tmp[i * 3 + 2];
		//printf("        Component: %d    H: %d    V: %d    Tq: %d\n", sof.nf[i].c, sof.nf[i].h, sof.nf[i].v, sof.nf[i].tq);
	}
}
void huffman(int id){
	uint16_t init = 0;
	int tmp = 0, count = 0;
	for(int i = 0; i < 16; i++){
		for(int j = 0; j < dht[id].li[i]; j++){
			dht[id].codeword[count] = init++;
			dht[id].codelength[count++] = i + 1;
			//printf("        (%2d,    %16d,       0x%02x)\n", dht[id].codelength[count - 1], dht[id].codeword[count - 1], dht[id].vi[count - 1]);
		}
		init <<= 1;
	}
}
void parseDHT(int id){
	int length;
	//printf("%x%x - DHT\n", buffer[0], buffer[1]);
	if(fread(buffer, 1, 2, fptr)){
		dht[id].length = buffer[0]*256 + buffer[1];
		length = dht[id].length - 2;
		//printf("        Lh: %d\n", dht[id].length);
	}
	while(id < 4){
		if(fread(buffer, 1, 1, fptr)){
			dht[id].tc = buffer[0] >> 4, dht[id].th = buffer[0] % 16;
			length--;
			//printf("        Tc: %d   Th: %d\n        Li: ", dht[id].tc, dht[id].th);
		}
		fread(dht[id].li, 16, 1, fptr); 
		int count = 0;
		for(int i = 0; i < 16; i++){
			//printf("%d ", dht[id].li[i]);
			count += dht[id].li[i];
			dht[id].accu[i] = count;
		}
		dht[id].size = count;
		//printf("\n        Vi,j: ");
		dht[id].vi = (uint8_t *)malloc(sizeof(uint8_t) * count);
		dht[id].codeword = (uint16_t *)malloc(sizeof(uint16_t) * count);
		dht[id].codelength = (uint16_t *)malloc(sizeof(uint16_t) * count);
		fread(dht[id].vi, count, 1, fptr); 
		for(int i = 0; i < count; i++)
			//printf("%d ", dht[id].vi[i]);
		//printf("\n        Size    CodeWord               Symbol\n");
		huffman(id);
		if(length - count - 16 > 0){
			length = length - count - 16, id++;
		}
		else 
			break;
	}
	
}
void print_bin(uint64_t a){
	for(int i = 0; i < 64; i++){
		if(i % 4 == 0)
			printf(" ");
		printf("%d", (a << i >> 63));
	}
	printf("\n");
}
void combine(uint64_t buffer64[], int bitcount, int* rest){
	int ret = 0;
	if(*rest <= bitcount){
		buffer64[0] = ((buffer64[0] << bitcount) + (buffer64[1] >> (64 - bitcount)));
		uint8_t tmp[8];
		for(int i = 0; i < 8; i++){
			fread(&tmp[i], 1, 1, fptr);
			if(tmp[i] == 0xFF){
				uint8_t gar;
				fread(&gar, 1, 1, fptr);
			}
			buffer64[1] = (buffer64[1] << 8) + tmp[i];
		}
		if(*rest != bitcount){
			buffer64[0] += (buffer64[1] >> (64 - (bitcount - *rest)));
			buffer64[1] <<= (bitcount - *rest);
		}
		*rest = 64 - (bitcount - *rest);
	}
	else{
		buffer64[0] = (buffer64[0] << bitcount) + (buffer64[1] >> (64 - bitcount));
		buffer64[1] <<= bitcount;
		*rest -= bitcount;	
	}
}
int extend(int T, int diff){
	if(diff < (1 << (T - 1)))
		return (1 - (1 << T) + diff);
	else
		return diff;
}
void IDCT(int* dest, int* source){
	float pi = 3.14159265;
	for(int i = 0; i < 8; i++){
		for(int j = 0; j < 8; j++){
			float sum = 0;
			for(int k = 0; k < 8; k++){
				for(int l = 0; l < 8; l++){
					//sum += cu[k] * cv[l] * cos(((2*i+1)*k*pi)/16) * cos(((2*j+1)*l*pi)/16) * source[k*8+l];	
					sum += cu[k] * cv[l] * cos_tb[i][k] * cos_tb[j][l] * source[k*8+l];	
				}
			}
			dest[i*8+j]  = sum/4;
		}
	}
}
int R(int y, int cb, int cr){
	float sum = y+1.402*(cr-128);
	if(sum >= 255)
		return 255;
	else if(sum <= 0)
		return 0;
	else
		return (uint8_t)sum;
}
int G(int y, int cb, int cr){
	float sum = y-0.34414*(cb-128)-0.71414*(cr-128);
	if(sum >= 255)
		return 255;
	else if(sum <= 0)
		return 0;
	else
		return (uint8_t)sum;
}
int B(int y, int cb, int cr){
	float sum = y+1.772*(cb-128);
	if(sum >= 255)
		return 255;
	else if(sum <= 0)
		return 0;
	else
		return (uint8_t)sum;
}
void BMP(YCbCr* ycbcr, int y, int x){
	int yh = sof.nf[0].h, yv = sof.nf[0].v;
	int cbh = sof.nf[1].h, cbv = sof.nf[1].v;
	int crh = sof.nf[2].h, crv = sof.nf[2].v;
	for(int i = 0; i < yv; i++){
		for(int j = 0; j < yh; j++){
			int yID = i*yh+j, ycbh = yh/cbh, ycbv = yv/cbv, ycrh = yh/crh, ycrv = yv/crv;
			int cbmax = ycbh*ycbv, crmax = ycrh*ycrv;
			int cbID = yID/cbmax, crID = yID/crmax;
			for(int k = 0; k < 8; k++){
				for(int l = 0; l < 8; l++){
					int cbi = (i*8+k)/ycbv, cbj = (j*8+l)/ycbh, cri = (i*8+k)/ycrv, crj = (j*8+l)/ycrh;
//					//printf("cbi = %d, cbj = %d------------------------\n", cbi, cbj);
					bmp[(y*yv+i)*8+k][(x*yh+j)*8+l][2] = R(ycbcr[0].idct[yID][k][l], ycbcr[1].idct[cbID][cbi][cbj],ycbcr[2].idct[crID][cri][crj]);
					bmp[(y*yv+i)*8+k][(x*yh+j)*8+l][1] = G(ycbcr[0].idct[yID][k][l], ycbcr[1].idct[cbID][cbi][cbj],ycbcr[2].idct[crID][cri][crj]);
					bmp[(y*yv+i)*8+k][(x*yh+j)*8+l][0] = B(ycbcr[0].idct[yID][k][l], ycbcr[1].idct[cbID][cbi][cbj],ycbcr[2].idct[crID][cri][crj]);
		//			//printf("Y %d, cb %d, cr %d--------------------------------\n", ycbcr[0].idct[yID][k][l], ycbcr[1].idct[cbID][cbi][cbj],ycbcr[2].idct[crID][cri][crj]);
//					//printf("R = %d, G = %d, B = %d\n",bmp[(y*yh+i)*8+k][(x*yv+j)*8+l][0], bmp[(y*yh+i)*8+k][(x*yv+j)*8+l][1], bmp[(y*yh+i)*8+k][(x*yv+j)*8+l][2]);
				}
			}
		}
	}
}
void decode(){
	uint8_t tmp[16];
	for(int i = 0; i < 16; i++){
		fread(&tmp[i], 1, 1, fptr);
		if(tmp[i] == 0xFF){
			uint8_t garbage;
			fread(&garbage, 1, 1, fptr);
		}
	}
	uint64_t buffer64[2];
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 8; j++){
			buffer64[i] = (buffer64[i] << 8) + tmp[i * 8 + j];
		}
	}
	int rest = 64, yDCPredict = 0, cbDCPredict = 0, crDCPredict = 0, predict[3];
	predict[0] = predict[1] = predict[2] = 0;
	maxY = sof.y / (8 * sof.nf[0].v), maxX = sof.x / (8 * sof.nf[0].h);
	if(sof.y % (8*sof.nf[0].v) != 0)
		maxY++;
	if(sof.x % (8*sof.nf[0].h) != 0)
		maxX++;
	for(int i = 0; i < maxY; i++){
		for(int j = 0; j < maxX; j++){
			//printf("MCU (%d, %d)\n", i, j);
			YCbCr ycbcr[3];
			memset((void *)ycbcr, 0, sizeof(YCbCr)*3);
			for(int k = 0; k < sof.nfSize; k++){
				for(int l = 0; l < sof.nf[k].v * sof.nf[k].h; l++){
					//printf("        DataUnit: (%d, %d, %d)\n", k, l >> 1, l % 2);
					
					//printf("                 DC Predictor:  %d\n", predict[k]);
					//printf("                 DC:\n");
					int DC = 0;
					int bit, id;
					for (int m = 0; m < 4; m++){
						if(dht[m].tc == 0 && sos.ns[k].td == dht[m].th)
							id = m;
					}
					int index = 0;
					int bitcount = dht[id].codelength[index];
					while(!DC){
						bit = buffer64[0] >> (64-bitcount);
						if(bit == dht[id].codeword[index]){
							combine(buffer64, bitcount, &rest);
							bitcount = dht[id].vi[index];
							if(bitcount == 0)
								bit = 0;
							else{
								bit = buffer64[0] >> (64 - bitcount);
								combine(buffer64, bitcount, &rest);
							}
							int ext = extend(dht[id].vi[index], bit);
							predict[k] += ext;
							//printf("                          T: %2d    DIFF: %2d       EXTEND(DIFF, T): %3d\n",dht[id].vi[index], bit, ext);
							ycbcr[k].table[l][0] = predict[k];
							DC = 1;
						}
						else{
							bitcount += dht[id].codelength[index+1] - dht[id].codelength[index];
							index++;
						}
					}
					
					//printf("                 AC:\n");
					for (int m = 0; m < 4; m++){
						if(dht[m].tc == 1 && sos.ns[k].ta == dht[m].th)
							id = m;
					}
					int EOB = 0, max = 1;
					while(!EOB){
						index = 0;
						bitcount = dht[id].codelength[index];
						int AC = 0;
						while(!AC){
							bit = buffer64[0] >> (64-bitcount);
							if(bit == dht[id].codeword[index]){
								combine(buffer64, bitcount, &rest);
								if(dht[id].vi[index] == 0){
									EOB = 1;
									break;
								}
								int run = dht[id].vi[index] >> 4, length = dht[id].vi[index] % 16;
								max += run;
								bitcount = length;
								if(bitcount == 0)
									bit = 0;
								else{
									bit = buffer64[0] >> (64 - bitcount);
									combine(buffer64, bitcount, &rest);
								}
								int ext= extend(length, bit);
								//printf("                          RS: 0x%02x RRRR: %2d SSSS: %2d ZZ(K): %2d EXTEND(ZZ(K), SSSS): %3d\n", dht[id].vi[index], run, length, bit, ext);
								ycbcr[k].table[l][max++] = ext;
								if(max >= 64){
									EOB = 1;
									break;
								}
								AC = 1;
							}
							else{
								bitcount += dht[id].codelength[index+1] - dht[id].codelength[index];
								index++;
							}
						}
					}
					zigzag(&ycbcr[k].zz[l][0][0], ycbcr[k].table[l], max);
					//printf("                          RS: 0x00 -- EOB\n"); 
				//	//printf("                 *Before Dequantize*\n"); 
					for(int m = 0; m < 8; m++){
				//		//printf("\n                 ");
						for(int n = 0; n < 8; n++){
				//			//printf("%4d ", ycbcr[k].zz[l][m][n]);
						}
					}
				//	//printf("\n");
				//	//printf("                 *Before IDCT*\n");
					int dqtID = (k > 1) ? 1 : k;
					for (int m = 0; m < 4; m++){
						if(sof.nf[k].tq == dqt[m].tq){
							dqtID = m;
							break;
						}
					}
					for(int m = 0; m < 8; m++){
				//		//printf("\n                 ");
						for(int n = 0; n < 8; n++){
							ycbcr[k].zz[l][m][n] *= dqt[dqtID].zz[m*8+n];
				//			//printf("%4d ", ycbcr[k].zz[l][m][n]);
						}
					}
				//	//printf("\n");
					IDCT(&ycbcr[k].idct[l][0][0], &ycbcr[k].zz[l][0][0]);
				//	//printf("                 *After IDCT*\n"); 
					for(int m = 0; m < 8; m++){
				//		//printf("\n                 ");
						for(int n = 0; n < 8; n++){
				//			//printf("%4d ", ycbcr[k].idct[l][m][n]);
						}
					}
				//	//printf("\n");
				//	//printf("                 *After +128*\n"); 
					for(int m = 0; m < 8; m++){
				//		//printf("\n                 ");
						for(int n = 0; n < 8; n++){
							ycbcr[k].idct[l][m][n] += (1<<7);
				//			//printf("%4d ", ycbcr[k].idct[l][m][n]);
						}
					}
				//	//printf("\n");

				}
			}
			BMP(ycbcr, i, j);
		}
	}
	//printf("0xffd9 - EOI\n");
	output(outputS);
}
void parseSOS(){
	//printf("%x%x - SOS\n", buffer[0], buffer[1]);
	if(fread(buffer, 1, 2, fptr)){
		sos.length = buffer[0]*256 + buffer[1];
		//printf("        Ls: %d   ", sos.length);
	}
	if(fread(buffer, 1, 1, fptr)){
		sos.nsSize = buffer[0];
		//printf("Ns: %d\n", sos.nsSize);
	}
	sos.ns = (NS *)malloc(sizeof(NS) * sos.nsSize);
	uint8_t tmp[2 * sos.nsSize];
	fread(tmp, 2 * sos.nsSize, 1, fptr);
	for(int i = 0; i < sos.nsSize; i++){
		sos.ns[i].cs = tmp[i * 2], sos.ns[i].td = tmp[i * 2 + 1] >> 4, sos.ns[i].ta = tmp[i * 2 + 1] % 16;
		//printf("        Cs: %d    Td: %d    Ta: %d\n", sos.ns[i].cs, sos.ns[i].td, sos.ns[i].ta);
	}
	fread(tmp, 3, 1, fptr);
	sos.ss = tmp[0], sos.se = tmp[1], sos.ah = tmp[2] >> 4, sos.al = tmp[2] % 16;
		//printf("        Ss: %d    Se: %d   Ah: %d    Al: %d\n", sos.ss, sos.se, sos.ah, sos.al);
	decode();	
}
void skip(){
	//printf("%x%x - skip\n", buffer[0], buffer[1]);
	if(fread(buffer, 1, 2, fptr)){
		int length = buffer[0]*256 + buffer[1];
		fseek ( fptr , length - 2 , SEEK_CUR );
	}
}
void parse(char * filename){
	fptr = fopen(filename, "rb");
	strcpy(outputS, filename);
	int lengthS = strlen(outputS);
	outputS[lengthS-3] = 'b', outputS[lengthS-2] = 'm', outputS[lengthS-1] = 'p';
	//printf("BMP\n");
	if(fptr != NULL){
		if(fread(buffer, 1, 2, fptr))
			//printf("%x%x - SOF\n", buffer[0], buffer[1]);
		fread(buffer, 1, 2, fptr); 
		while(buffer[0]*256 + buffer[1] != EOI_SYMB){
			switch(buffer[0]*256 + buffer[1]){
				case APP_SYMB:
					parseAPP();
					break;
				case DQT_SYMB:
					parseDQT(dqt_count++);
					break;
				case SOF_SYMB:
					parseSOF();
					break;
				case DHT_SYMB:
					parseDHT(dht_count++);
					break;
				case SOS_SYMB:
					parseSOS();
					break;
				default:
					skip();
					break;
			}
			fread(buffer, 1, 2, fptr);
		}
	}
	fclose(fptr);
}
void output(char * filename){
	typedef struct { 
		/* type : Magic identifier,一般為BM(0x42,0x4d) */ 
		unsigned short int type; 
		unsigned int size;/* File size in bytes,全部的檔案大小 */ 
		unsigned short int reserved1, reserved2; /* 保留欄位 */ 
		unsigned int offset;/* Offset to image data, bytes */ 
	} FILEHEADER;
	typedef struct { 
		unsigned int size;/* Info Header size in bytes */ 
		int width,height;/* Width and height of image */ 
		unsigned short int planes;/* Number of colour planes */ 
		unsigned short int bits; /* Bits per pixel */ 
		unsigned int compression; /* Compression type */ 
		unsigned int imagesize; /* Image size in bytes */ 
		int xresolution,yresolution; /* Pixels per meter */ 
		unsigned int ncolours; /* Number of colours */ 
		unsigned int importantcolours; /* Important colours */ 
	} INFOHEADER;
	FILEHEADER fileheader;
	int Y = sof.y, X = (sof.x*3+3)/4*4;
	//int Y = maxY*sof.nf[0].v*8, X = maxX*sof.nf[0].h*8;
	fileheader.type = 0x4d42, fileheader.size = 54+X*Y, fileheader.reserved1 = fileheader.reserved2 = 0, fileheader.offset = 54;
	INFOHEADER infoheader;
	infoheader.size = 40, infoheader.width = sof.x, infoheader.height = sof.y, infoheader.planes = 1, infoheader.bits = 24, infoheader.compression = 0;
	infoheader.imagesize = 0, infoheader.xresolution = infoheader.yresolution = 0x0b12, infoheader.ncolours = infoheader.importantcolours = 0;	
	FILE * wfp;
	wfp = fopen(filename, "wb");
	fwrite(&fileheader.type, 1, sizeof(fileheader.type), wfp);
	fwrite(&fileheader.size, 1, sizeof(fileheader.size), wfp);
	fwrite(&fileheader.reserved1, 1, sizeof(fileheader.reserved1), wfp);
	fwrite(&fileheader.reserved2, 1, sizeof(fileheader.reserved2), wfp);
	fwrite(&fileheader.offset, 1, sizeof(fileheader.offset), wfp);
	fwrite(&infoheader, 1, sizeof(INFOHEADER), wfp);
	uint8_t imgbuf[X];
	for(int i = Y-1; i >= 0; i--){
		memset(imgbuf, 0, sizeof(imgbuf));
		for(int j = 0; j < sof.x; j++){
			imgbuf[j*3] = bmp[i][j][0];
			imgbuf[j*3+1] = bmp[i][j][1];
			imgbuf[j*3+2] = bmp[i][j][2];
		}
		fwrite(imgbuf, sizeof(uint8_t), X, wfp);
	}
	
	fclose(wfp);	
	exit(0);
}
