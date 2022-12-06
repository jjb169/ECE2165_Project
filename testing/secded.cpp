//james bickerstaff
//testing basic double data type functions for SEC-DED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>

#define DEBUG 0

/* P matrix and syndromes created by https://github.com/jgaeddert/liquid-dsp/blob/master/src/fec/src/fec_secded7264.c*/

// P matrix [8 x 64] - H matrix if you add on the Identity matrix
//  11111111 00001111 00001111 00001100 01101000 10001000 10001000 10000000 : 1 0 0 0 0 0 0 0 ; t65=s63^...^s56^s51^...^s48^s43^...^s40^s35^s34^s30^s29^s27^s23^s19^s15^s11^s7
//  11110000 11111111 00000000 11110011 01100100 01000100 01000100 01000000 : 0 1 0 0 0 0 0 0 ; t66=s63^...^s60^s55...^s48^s39^...^s36^s33^s32^s30^s29^s26^s22^s18^s14^s10^s6
//  00110000 11110000 11111111 00001111 00000010 00100010 00100010 00100110 : 0 0 1 0 0 0 0 0 ; t67=s61^s60^s55^...^s52^s47^...^s40^s35^...^s32^s25^s21^s17^s13^s9^s5^s2^s1
//  11001111 00000000 11110000 11111111 00000001 00010001 00010001 00010110 : 0 0 0 1 0 0 0 0 ; t68=s63^s62^s59^...^56^s47^...^s44^s39^...^s32^s24^s20^s16^s12^s8^s4^s2^s1
//  01101000 10001000 10001000 10000000 11111111 00001111 00000000 11110011 : 0 0 0 0 1 0 0 0 ; t69=s62^s61^s59^s55^s51^s47^s43^s39^s31^...^s24^s19^...^s16^s7^...^s4^s2^s1
//  01100100 01000100 01000100 01000000 11110000 11111111 00001111 00001100 : 0 0 0 0 0 1 0 0 ; t70=s62^s61^s58^s54^s50^s46^s42^s38^s31^...^s28^s23^...^s16^s11^...^s8^s3^s2
//  00000010 00100010 00100010 00100110 11001111 00000000 11111111 00001111 : 0 0 0 0 0 0 1 0 ; t71=s57^s53^s49^s45^s41^s37^s34^s33^s31^s30^s27^...^s24^s15^...^s8^s3^...^s0
//  00000001 00010001 00010001 00010110 00110000 11110000 11110000 11111111 : 0 0 0 0 0 0 0 1 ; t72=s56^s52^s48^s44^s40^s36^s34^s33^s29^s28^s23^...^s20^s15^...^s12^s7^...^s0
unsigned char P_matrix[64] = {
    0xFF, 0x0F, 0x0F, 0x0C, 0x68, 0x88, 0x88, 0x80,
    0xF0, 0xFF, 0x00, 0xF3, 0x64, 0x44, 0x44, 0x40,
    0x30, 0xF0, 0xFF, 0x0F, 0x02, 0x22, 0x22, 0x26,
    0xCF, 0x00, 0xF0, 0xFF, 0x01, 0x11, 0x11, 0x16,
    0x68, 0x88, 0x88, 0x80, 0xFF, 0x0F, 0x00, 0xF3,
    0x64, 0x44, 0x44, 0x40, 0xF0, 0xFF, 0x0F, 0x0C,
    0x02, 0x22, 0x22, 0x26, 0xCF, 0x00, 0xFF, 0x0F,
    0x01, 0x11, 0x11, 0x16, 0x30, 0xF0, 0xF0, 0xFF};

// syndrome vectors for errors of weight 1 (1 error exists) - each index is a syndrome
unsigned char syndromes[72] = {
    0x0b, 0x3b, 0x37, 0x07, 0x19, 0x29, 0x49, 0x89,
    0x16, 0x26, 0x46, 0x86, 0x13, 0x23, 0x43, 0x83,
    0x1c, 0x2c, 0x4c, 0x8c, 0x15, 0x25, 0x45, 0x85,
    0x1a, 0x2a, 0x4a, 0x8a, 0x0d, 0xcd, 0xce, 0x0e,
    0x70, 0x73, 0xb3, 0xb0, 0x51, 0x52, 0x54, 0x58,
    0xa1, 0xa2, 0xa4, 0xa8, 0x31, 0x32, 0x34, 0x38,
    0xc1, 0xc2, 0xc4, 0xc8, 0x61, 0x62, 0x64, 0x68,
    0x91, 0x92, 0x94, 0x98, 0xe0, 0xec, 0xdc, 0xd0,
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
	

//method to append check byte/char
void appendCheck(char data[])
{
	//must do the parity checks of the very long XOR strings above against the data string
	char checkBits[8];
	//parity variable, set as 0 first due to even parity
	char parity = 0;
	//this for loop to assign checkBits 0 - 7
	for(int i = 0; i < 8; i++)
	{
		//this for loop to go through a row of the P matrix and 8 data bytes
		for(int d_index = 8, p_index = 0; d_index > 1; d_index--, p_index++)
		{
			//grab the current data byte and mask it with the P index to check this byte's parity
			char workingVar = data[d_index] & P_matrix[p_index+(8*i)];
			
			if(DEBUG == 1)
			{
				//extreme debugging
				printf("\nWorking var = %02X\n", workingVar);
				for(int j = 7; j >= 0; j--)
				{
					//std::cout << "curr = " << (*(bytes+i) >> j) << "\n";
					if(workingVar >> j & 1 == 1)
						printf("1");
					else
						printf("0");
				}
			}
			
			
			//loop through this byte and count 1s
			for(int j = 0; j < 8; j++)
			{
				//check if current bit is a 1
				if((workingVar >> j) & 1 == 1)
				{
					//if yes, increment parity counter mod 2
					parity = (parity + 1) % 2;
				}
			}
			
			//checkBits[i] = checkBits[i] ^ (data[d_index] ^ P_matrix[p_index+(8*i)])
		}
		
		if(DEBUG == 1)
		{
			printf("\nparity = %d", parity);
			printf("\nNext check bit\n");
		}
		//assign this parity to the check bit
		checkBits[i] = parity;
		//reset parity for next checkBit 
		parity = 0;
	}
	
	if(DEBUG == 1)
	{
		//print check bits
		printf("\nCheck bits = "); 
		for(int i = 0; i < 8; i++)
		{
			printf("%01X ", checkBits[i]); 
		}
		printf("\n");
	}

	//final check byte
	char checkByte = checkBits[0] | (checkBits[1] << 1) | (checkBits[2] << 2) | (checkBits[3] << 3) | (checkBits[4] << 4) | (checkBits[5] << 5) | (checkBits[6] << 6) | (checkBits[7] << 7);
	
	//append check bit onto data bit to create whole package
	data[0] = checkByte;
}

//method to produce the syndrome from received data
void calculateSyndrome(char data[])
{
	//must do the parity checks of the very long XOR strings above against the data string
	char checkBits[8];
	//parity variable, set as 0 first due to even parity
	char parity = 0;
	//this for loop to assign checkBits 0 - 7
	for(int i = 0; i < 8; i++)
	{
		//this for loop to go through a row of the P matrix and 8 data bytes
		for(int d_index = 8, p_index = 0; d_index > 1; d_index--, p_index++)
		{
			//grab the current data byte and mask it with the P index to check this byte's parity
			char workingVar = data[d_index] & P_matrix[p_index+(8*i)];
			
			if(DEBUG == 1)
			{
				//extreme debugging
				printf("\nWorking var = %02X\n", workingVar);
				for(int j = 7; j >= 0; j--)
				{
					//std::cout << "curr = " << (*(bytes+i) >> j) << "\n";
					if(workingVar >> j & 1 == 1)
						printf("1");
					else
						printf("0");
				}
			}
			
			
			//loop through this byte and count 1s
			for(int j = 0; j < 8; j++)
			{
				//check if current bit is a 1
				if((workingVar >> j) & 1 == 1)
				{
					//if yes, increment parity counter mod 2
					parity = (parity + 1) % 2;
				}
			}
			
			//checkBits[i] = checkBits[i] ^ (data[d_index] ^ P_matrix[p_index+(8*i)])
		}
		
		//finally, check if this lines parity/debug bit is a 1
		if((data[0] >> i) & 1 == 1)
		{
			//if yes, increment parity counter mod 2
			parity = (parity + 1) % 2;
		}
		
		if(DEBUG == 1)
		{
			printf("\nparity = %d", parity);
			printf("\nNext check bit\n");
		}
		//assign this parity to the check bit
		checkBits[i] = parity;
		//reset parity for next checkBit 
		parity = 0;
	}
	
	if(DEBUG == 1)
	{
		//print check bits
		printf("\nCheck bits = "); 
		for(int i = 0; i < 8; i++)
		{
			printf("%01X ", checkBits[i]); 
		}
		printf("\n");
	}

	//attach all bits to create the syndrome
	char syndrome = checkBits[0] | (checkBits[1] << 1) | (checkBits[2] << 2) | (checkBits[3] << 3) | (checkBits[4] << 4) | (checkBits[5] << 5) | (checkBits[6] << 6) | (checkBits[7] << 7);
	printf("Syndrome produced is: %02X\n", syndrome);
	
}



//main function to test external methods
int main(int argc, char** argv)
{
	//take a double as input from keyboard for testing
	//actual variable
	double test;
	char *bytes;
	printf("Please input a number to test: ");
	std::cin >> test;
	bytes = (char *)(&test);	
	
	//print number as output
	printf("The number was: %.3lf\n", test);
	
	//print the number in binary
	printf("The number in binary/hex is: ");
	for(int i = 7; i >= 0; i--)
	{
		printf("\nByte %d = %02X\n", i, *(bytes+i));
		for(int j = 7; j >= 0; j--)
		{
			//std::cout << "curr = " << (*(bytes+i) >> j) << "\n";
			if(*(bytes+i) >> j & 1 == 1)
				printf("1");
			else
				printf("0");
		}
	}
	printf("\n\n");
	
	//message to transmit: will be a char array[9] where [0] is check bits and [1-8] are data bits
	char transmit[9];
	transmit[0] = 0;
	//assign data bytes now
	for(int i = 1; i < 9; i++)
	{
		transmit[i] = *(bytes+i-1);
	}
	//print for debug
	printf("Printing transmit data before adding check bits: \n");
	for(int i = 8; i >= 0; i--)
	{
		printf("%02X",transmit[i]);
	}
	printf("\n\n");
	
	//need to determine the byte to add on through parity checks
	appendCheck(transmit);
	printf("Printing transmit data after adding check bits: \n");
	for(int i = 8; i >= 0; i--)
	{
		printf("%02X",transmit[i]);
	}
	printf("\n\n");
	
	
	//alter data in some way
	transmit[8] = transmit[8] | 0x80;
	printf("Printing transmit data after messing it up: \n");
	for(int i = 8; i >= 0; i--)
	{
		printf("%02X",transmit[i]);
	}
	printf("\n\n");
	
	//calculate syndrome
	calculateSyndrome(transmit);
	
	return 0;
}

