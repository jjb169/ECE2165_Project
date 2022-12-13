//james bickerstaff
//testing basic double data type functions for SEC-DED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <cstdlib>

#define DEBUG 0
#define BIT_DEBUG 0
#define PRINTINFO 0

/* P matrix and syndromes created by https://github.com/jgaeddert/liquid-dsp/blob/master/src/fec/src/fec_secded7264.c */

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
//   7     6     5     4     3     2     1     0

//method to append check byte/char
void appendCheck(unsigned char data[])
{
	//must do the parity checks of the very long XOR strings above against the data string
	char checkBits[8];
	//parity variable, set as 0 first due to even parity
	char parity = 0;
	//this for loop to assign checkBits 0 - 7
	for(int i = 0; i < 8; i++)
	{
		//this for loop to go through a row of the P matrix and 8 data bytes
		for(int d_index = 8, p_index = 0; d_index > 0; d_index--, p_index++)
		{
			//grab the current data byte and mask it with the P index to check this byte's parity
			char workingVar = data[d_index] & P_matrix[p_index+(8*i)];
			
			if(BIT_DEBUG == 1)
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
		
		if(BIT_DEBUG == 1)
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
char calculateSyndrome(unsigned char data[])
{	
	//must do the parity checks of the very long XOR strings above against the data string
	char checkBits[8];
	//parity variable, set as 0 first due to even parity
	char parity = 0;
	//this for loop to assign checkBits 0 - 7
	for(int i = 0; i < 8; i++)
	{
		//this for loop to go through a row of the P matrix and 8 data bytes
		for(int d_index = 8, p_index = 0; d_index > 0; d_index--, p_index++)
		{
			//grab the current data byte and mask it with the P index to check this byte's parity
			char workingVar = data[d_index] & P_matrix[p_index+(8*i)];
			
			if(BIT_DEBUG == 1)
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
		
		if(BIT_DEBUG == 1)
		{
			printf("\nparity = %d", parity);
			printf("\nNext check bit\n");
		}
		//assign this parity to the check bit
		checkBits[i] = parity;
		//reset parity for next checkBit 
		parity = 0;
	}
	
	if(BIT_DEBUG == 1)
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
	unsigned char syndrome = checkBits[0] | (checkBits[1] << 1) | (checkBits[2] << 2) | (checkBits[3] << 3) | (checkBits[4] << 4) | (checkBits[5] << 5) | (checkBits[6] << 6) | (checkBits[7] << 7);
	if(PRINTINFO)
	{
		printf("Syndrome produced is: %02X\n", syndrome);
	}
	//printf("Syndrome produced unmasked is: %02X\n", syndrome );
	return syndrome;
}


//method to find the incorrect bit and correct if possible
bool secded_data(unsigned char data[], unsigned char syndrome)
{
	//ints to store result before fixing
	int byteError = 8;
	int bitError = 0;
	int bitNum = 0;
	//bool for checking if a matching syndrome is found
	bool syndromeMatch = false;
	//printf("Trying to find: %02X\n", syndrome);
	//first check syndrome against list of syndromes to find bit in error
	for(int i = 71; i >= 0; i--)
	{
		//printf("Syndromes[i] = %02X\n", syndromes[i]);
		//check curr syndrome for match
		if(syndrome == syndromes[i])
		{
			//printf("Syndrome match: %02X\n", syndromes[i]);
			//printf("i = %d, bitNum = %d\n", i, bitNum); 
			byteError = bitNum / 8;
			bitError = bitNum % 8;
			//we now know which byte and bit are in error
			syndromeMatch = true;
			break;
		}
		//since syndrome matrix is stored in reverse order, need a separate variable for when doing correction
		bitNum++;
	}
	
	//I can't seem to figure out the issue with correcting the check bits themselves, not a big deal and prob dont actually need corrected
	//but I'll do it anyways
	//appears that the exact opposite bit is identified as being in error, ie 0 when it should be 7, 4 when it should be 3, etc, so flip em
	//very hacky
	if(byteError == 0)
	{
		if(bitError == 0)
			bitError = 7;
		else if(bitError == 1)
			bitError = 6;
		else if(bitError == 2)
			bitError = 5;
		else if(bitError == 3)
			bitError = 4;
		else if(bitError == 4)
			bitError = 3;
		else if(bitError == 5)
			bitError = 2;
		else if(bitError == 6)
			bitError = 1;
		else if(bitError == 7)
			bitError = 0;
	}
	
	if(DEBUG == 1)
	{
		printf("Byte in error is: %d\n", byteError);
		printf("Bit in error is: %d\n", bitError);
	}
	
	//correct error if possible
	if(syndromeMatch)
	{
		//if 0 flip to 1
		if(data[byteError] & (0x01 << bitError) == 0)
		{
			//printf("Correcting 0 to 1\n");
			data[byteError] = data[byteError] | (0x01 << bitError);
		}
		else
		{
			//else it is a 0, flip to 1
			//printf("Correcting 1 to 0\n");
			data[byteError] = data[byteError] ^ (0x01 << bitError);
		}	
	}
	
	//return syndrome match to signify if bit has been corrected or if double error is corrected
	return syndromeMatch;
}


//function to take in a double, add the check byte, and place the results in a passed in char array
//inputs:
// 1) a double that contains the data being transmitted
// 2) an UNSIGNED char array of NINE (9) indices to hold 8 bytes of data and 1 byte of check (64 bits + 8 = 72 bits for (72,64) Hamming code)
//output:
// just modifying the input char array directly, no return values
void addCheckByte(double data, unsigned char codeWord[])
{	
	//turn the input data into a byte array
	unsigned char *dataBytes = (unsigned char *)(&data);	
	
	//message to transmit: will be a char array[9] where [0] is check bits and [1-8] are data bits
	codeWord[0] = 0;
	//assign data bytes now
	for(int i = 1; i < 9; i++)
	{
		codeWord[i] = *(dataBytes+i-1);
	}
	
	//need to determine the byte to add on through parity checks
	appendCheck(codeWord);
	
	//parity byte added to char array, ready to transmit
}


//function to modify bits - currently only 1 or 2
//inputs:
// 1) codeWord in the form of unsigned char array -> 9 bytes
// 2) chance of a SINGLE bit error occurring, just using an int for "X mod chance" calculation
// 3) chance of a DOUBLE bit error occurring, just using an int for "X mod chance" calculation
//output:
// returning the number of modified bits
// modifying the codeWord array directly to induce bit errors
int modifyBits(unsigned char codeWord[], int singleErrorChance, int doubleErrorChance)
{
	//var to return how many bits were modified
	int modifiedCount = 0; 
	
	//alter data in some way
	//check if random value generated mod error chance is 0, if so, add a single bit error
	if(std::rand() % singleErrorChance == 0)
	{
		//increment bit error count
		modifiedCount++;
		//first select which byte and bit to modify
		int byte = std::rand() % 9; //byte to alter
		int bit = std::rand() % 8; //bit withihn that byte
		//if 0, change to 1, else it is a 1, change to 0
		if(codeWord[byte] & (0x01 << bit) == 0)
		{
			codeWord[byte] = codeWord[byte] | (0x01 << bit);
		}
		else
		{
			codeWord[byte] = codeWord[byte] ^ (0x01 << bit);
		}
		
		//determine if doing double error as well, doing mod doubleErrorChance 
		if(std::rand() % doubleErrorChance == 0)
		{
			//increment bit error count
			modifiedCount++;
			
			//save old byte/bit values so we bit flip something else instead of flipping it back to normal
			int byteTemp = byte, bitTemp = bit;
			
			//once again generate a byte and bit to mess with
			while(byteTemp == byte)
				byte = std::rand() % 9; //byte to alter
			while(bitTemp == bit)
				bit = std::rand() % 8; //bit withihn that byte
			//if 0, change to 1, else it is a 1, change to 0
			if(codeWord[byte] & (0x01 << bit) == 0)
			{
				codeWord[byte] = codeWord[byte] | (0x01 << bit);
			}
			else
			{
				codeWord[byte] = codeWord[byte] ^ (0x01 << bit);
			}
		}
	}
	
	//return # of modified bits
	return modifiedCount;	
}


//modify to check - and possibly correct the received data after a transmission
//input:
// 1) receivedData is the data from a transmission that may or may not be modified, in form of unsigned char array
//output:
// bool signifying whether the data is correct or could not be fixed
// also directly fix the data array of received data
bool checkReceivedData(unsigned char receivedData[])
{	
	//bool for if data has been fixed (or not needed fixed) or needs retransmitted (DED)
	bool correctData = false;
	//calculate syndrome based on received data
	unsigned char syndrome = calculateSyndrome(receivedData);
	
	//do all the following if syndrome != 0, meaning no error
	if(syndrome != 0x00)
	{
		//correct bit if possible
		correctData = secded_data(receivedData, syndrome);
		
		//if it could not be corrected,double error has been detected
		if(correctData)
		{
			//DATA SUCCESSFULLY FIXED :D
			//dont really need to do anything else, bool already set to true for return
			correctData = true;
		}
		else
		{
			//DED DETECTED, NEED TO RETRANSMIT DATA
			correctData = false;
		}		
	}
	else
	{
		//NO ERROR DETECTED, DATA NEEDS NO MODIFICATION
		correctData = true;
	}
	
	//return bool signifying if data is good or not
	return correctData;
}

//simple method to retranslate data from unsigned char array back to a double
//input:
// 1) unsigned char array containing 9 bytes of 8 data + 1 check byte
//output:
// returns the data received and fixed as a double, as needed by program
double charToDouble(unsigned char receivedData[])
{
	//want just bytes 1 through 8, as byte 0 is a check
	unsigned char originalData[8];
	for(int i = 1, j = 0; i < 9; i++, j++)
	{
		//copy byte by byte
		originalData[j] = receivedData[i];
	}
	//turn the original 8 bytes back into a double
	double fixedData = *(double*)originalData;
	
	//return this original data
	return fixedData;
}


//modified main function to make use of methods for each step, as will be done in actual program
int main(int argc, char** argv)
{
	//grab iteration count from cmd
	if(argc < 4 || argc > 4)
	{
		printf("Please have 3 arguments: # iterations for fault injection, mod for single error chance, mod for double error chance\n");
		exit(0);
	}
	//grab iteration count
	int totalIterations = std::atoi(argv[1]);
	int singleErrorChance = std::atoi(argv[2]);
	int doubleErrorChance = std::atoi(argv[3]);
	
	//seed rng
	std::srand(time(0));
	
	//PHASE 1: TAKE INPUT AND ADD CHECK BITS
	
	//take a double as input from keyboard for testing
	//actual variable
	double test;
	unsigned char *bytes;
	printf("Please input a number to test: ");
	std::cin >> test;
	bytes = (unsigned char *)(&test);	
	//message to transmit: will be a char array[9] where [0] is check bits and [1-8] are data bits
	unsigned char transmit[9];
	//send in double and char array to add check bits
	addCheckByte(test, transmit);
	
	//SETUP FOR LOOPING PHASES 2 AND 3 (FAULT INJECTION AND CORRECTION)
	
	//counting vars
	int se_corrections = 0, ded_count = 0, no_error_count = 0;
	int se_faults = 0, ded_faults = 0;
	int iters = 0;
	
	//do this X times for testing
	for(iters = 0; iters < totalIterations; iters++)
	{
		//copy data for this iteration
		unsigned char copy[9];
		for(int j = 0; j < 9; j++)
		{
			copy[j] = transmit[j];
		}
	
		//PHASE 2: INJECT A FAULT OR TWO
		
		int numBitsModified = modifyBits(copy, singleErrorChance, doubleErrorChance);
		//update counter variables for output verification
		if(numBitsModified == 1)
		{
			se_faults++;
		}
		else if(numBitsModified == 2)
		{
			ded_faults++;
		}
		
		//PHASE 3: TAKE IN MODIFIED DATA, CALCULATE SYNDROME AND DECIDE NEXT COURSE OF ACTION
		
		bool corrected = checkReceivedData(copy);
		
		//check result of correction attempt
		
		//4 cases:
		//corrected == true and numBitsModified == 0, meaning no error detected
		//corrected == true and numBitsModified == 1, meaning single bit error fixed
		//corrected == false and numBitsModified == 2, meaning double bit error detected, cannot correct
		//corrected == false and numBitsModified != 2, meaning we missed a correction
		if(corrected && numBitsModified == 0)
		{
			no_error_count++;
			//also want to turn data back into a double and make sure it matches original input
			double correctedData = charToDouble(copy);
			if(test != correctedData)
			{
				printf("NO MODIFICATIONS BUT TRANSLATED DATA DOES NOT MATCH ORIGINAL INPUT!\n");
			}
		}
		else if(corrected && numBitsModified == 1)
		{
			se_corrections++;
			//also want to turn data back into a double and make sure it matches original input
			double correctedData = charToDouble(copy);
			if(test != correctedData)
			{
				printf("CORRECTED DATA DOES NOT MATCH ORIGINAL INPUT!\n");
			}
		}
		else if (!corrected && numBitsModified == 2)
		{
			ded_count++;
		}
		else
		{
			//missed correcting or detecting an error
			printf("MISSED AN ERROR!!\n");
			printf("Original data: ");
			for(int i = 8; i >= 0; i--)
			{
				printf("%02X",transmit[i] & 0xFF);
			}
			printf("Modified data: ");
			for(int i = 8; i >= 0; i--)
			{
				printf("%02X",copy[i] & 0xFF);
			}
			printf("\n");
		}
	}
	
	//done testing, output results
	//print results
	printf("Out of %d iterations:\n", iters);
	printf("Single Errors injected: %d\n", se_faults);
	printf("Single Error corrections: %d\n", se_corrections);
	printf("Double Error detections: %d\n", ded_count);
	printf("Double Errors injected: %d\n", ded_faults);
	printf("No Errors found: %d\n", no_error_count);
}


/*
//main function to test external methods - DONT USE/ONLY USE FOR VERBOSE DEBUGGING PURPOSES
int mainBackup(int argc, char** argv)
{
	
	//grab iteration count from cmd
	if(argc < 4 || argc > 4)
	{
		printf("Please have 3 arguments: # iterations for fault injection, mod for single error chance, mod for double error chance\n");
		exit(0);
	}
	//grab iteration count
	int totalIterations = std::atoi(argv[1]);
	int singleErrorChance = std::atoi(argv[2]);
	int doubleErrorChance = std::atoi(argv[3]);
	
	//seed rng
	std::srand(time(0));
	
	
	//PHASE 1: TAKE INPUT AND ADD CHECK BITS ********************************************
	
	//take a double as input from keyboard for testing
	//actual variable
	double test;
	unsigned char *bytes;
	printf("Please input a number to test: ");
	std::cin >> test;
	bytes = (unsigned char *)(&test);	
	
	//print number as output
	printf("The number was: %.3lf\n", test);
	
	//print the number in binary
	printf("The number in binary/hex is: ");
	for(int i = 7; i >= 0; i--)
	{
		printf("\nByte %d = %02X\n", i, (*(bytes+i)) & 0xFF);
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
	unsigned char transmit[9];
	transmit[0] = 0;
	//assign data bytes now
	for(int i = 1; i < 9; i++)
	{
		transmit[i] = *(bytes+i-1);
	}
	//print for debug
	printf("Printing transmit data before modifying check bits (final 2 hex #s): \n");
	for(int i = 8; i >= 0; i--)
	{
		printf("%02X",transmit[i] & 0xFF);
	}
	printf("\n\n");
	
	//need to determine the byte to add on through parity checks
	appendCheck(transmit);
	printf("Printing transmit data after adding check bits (final 2 hex #s): \n");
	for(int i = 8; i >= 0; i--)
	{
		printf("%02X",transmit[i] & 0xFF);
	}
	printf("\n\n");
	
	//END PHASE 1
	
	//counting vars
	int se_corrections = 0, ded_count = 0, no_error_count = 0;
	int se_faults = 0, ded_faults = 0;
	int iters = 0;
	
	//do this X times for testing
	for(iters = 0; iters < totalIterations; iters++)
	{
		//PHASE 2: INJECT A FAULT OR TWO ********************************************************
		
		//bools for outputting missed single or double bit errors
		bool de_injected = false, se_injected = false;
		bool ded_missed = true, sec_missed = true;
		
		//copy data for this iteration
		unsigned char copy[9];
		for(int j = 0; j < 9; j++)
		{
			copy[j] = transmit[j];
		}
		//alter data in some way
		//single error likelihood, just doing mod 2 for ~50%
		if(std::rand() % singleErrorChance == 0)
		{
			//printf("Inserting a single bit error\n");
			//change var so we know a single error was injected
			se_injected = true;
			//increment single error fault count
			se_faults++;
			
			int byte = std::rand() % 9; //byte to alter
			int bit = std::rand() % 8; //bit withihn that byte
			//if 0, change to 1, else it is a 1, change to 0
			if(copy[byte] & (0x01 << bit) == 0)
			{
				if(PRINTINFO)
					printf("Changing 0 to 1\n");
				copy[byte] = copy[byte] | (0x01 << bit);
			}
			else
			{
				if(PRINTINFO)
					printf("Changing 1 to 0\n");
				copy[byte] = copy[byte] ^ (0x01 << bit);
			}
			
			//debug/printing information
			if(PRINTINFO)
			{
				printf("Single error:\n");
				printf("Changing byte %d, bit %d\n", byte, bit); 
				printf("Printing data after messing it up: \n");
				for(int i = 8; i >= 0; i--)
				{
					printf("%02X",copy[i] & 0xFF);
				}
				printf("\n");
			}
			
			//determine if doing double error as well, doing mod 10 for ~10% chance following the 50% chance of single error
			//so like 5% chance this occurs
			if(std::rand() % doubleErrorChance == 0)
			{
				//change var so we know a double error was injected
				de_injected = true;
				//increment double error fault count
				ded_faults++;
				
				//save old byte/bit values so we bit flip something else instead of flipping it back to normal
				int byteTemp = byte, bitTemp = bit;
				
				//once again generate a byte and bit to mess with
				while(byteTemp == byte)
					byte = std::rand() % 9; //byte to alter
				while(bitTemp == bit)
					bit = std::rand() % 8; //bit withihn that byte
				//if 0, change to 1, else it is a 1, change to 0
				if(copy[byte] & (0x01 << bit) == 0)
				{
					if(PRINTINFO)
						printf("Changing 0 to 1\n");
					copy[byte] = copy[byte] | (0x01 << bit);
				}
				else
				{
					if(PRINTINFO)
						printf("Changing 1 to 0\n");
					copy[byte] = copy[byte] ^ (0x01 << bit);
				}
				
				//printf("Inserting a double bit error\n");
				
				//debug/printing information
				if(PRINTINFO)
				{
					printf("Double error:\n");
					printf("Changing byte %d, bit %d\n", byte, bit); 
					printf("Printing data after messing it up: \n");
					for(int i = 8; i >= 0; i--)
					{
						printf("%02X",copy[i] & 0xFF);
					}
					printf("\n");
				}
			}
		}
		
		
		//PHASE 3: TAKE IN MODIFIED DATA, CALCULATE SYNDROME AND DECIDE NEXT COURSE OF ACTION ************
		
		//calculate syndrome
		unsigned char syndrome = calculateSyndrome(copy);
		
		//do all the following if syndrome != 0, meaning no error
		if(syndrome != 0x00)
		{
			//correct bit if possible
			bool corrected = secded_data(copy, syndrome);
			
			//if it could not be corrected,double error has been detected
			if(corrected)
			{
				//set sec_missed to false to indicate we caught this error
				sec_missed = false;
				
				//check if data matches that sent, verify if correction was done properly
				bool matches = true;
				for(int j = 0; j < 9; j++)
				{
					if(copy[j] != transmit[j])
					{
						matches = false;
					}
				}
				//increment success variable if done properly
				if(matches)
				{
					se_corrections++;
				}
				
				//print extra info
				if(PRINTINFO)
				{
					//printing data after fixing it and stating if it is a match with OG data
					printf("Data after correcting error: \n");
					for(int i = 8; i >= 0; i--)
					{
						printf("%02X",copy[i] & 0xFF);
					}
					printf("\n");
					//if it matches, state so, otherwise state error
					if(matches)
					{
						printf("Error corrected successfully! PASS\n\n");
					}
					else
					{
						printf("Error not corrected successfully! FAIL\n\n");
					}
				}
			}
			else
			{
				//make ded_missed var false to indicate we caught this error
				ded_missed = false;
				//increment ded variable
				ded_count++;
				//print extra output
				if(PRINTINFO)
				{
					//double error
					printf("Double error has been detected, not correctable!\n\n");
				}
			}
			
			
			//PRINT OUT DATA REGARDING MISSED FAULT - NOT NEEDED IN ACTUAL APP REALLY
			//display info regarding a missed fault
			if(se_injected)
			{
				//also check if it was a double error
				if(de_injected)
				{
					//check if we caught the double error
					if(ded_missed)
					{
						//display info regarding what was missed
						printf("MISSED A DOUBLE ERROR!!\n");
						printf("Original data: ");
						for(int i = 8; i >= 0; i--)
						{
							printf("%02X",transmit[i] & 0xFF);
						}
						printf("Modified data: ");
						for(int i = 8; i >= 0; i--)
						{
							printf("%02X",copy[i] & 0xFF);
						}
						printf("\n");
					}
				}
				else
				{
					//else, was a single error, display info about miss if applicable
					if(sec_missed)
					{
						printf("MISSED A SINGLE ERROR!!\n");
						printf("Original data: ");
						for(int i = 8; i >= 0; i--)
						{
							printf("%02X",transmit[i] & 0xFF);
						}
						printf("Modified data: ");
						for(int i = 8; i >= 0; i--)
						{
							printf("%02X",copy[i] & 0xFF);
						}
						printf("\n");
					}
				}
			}
			
		}
		else
		{
			//increment no error count
			no_error_count++;
			
			if(PRINTINFO)
			{
				//no error found
				printf("No error detected!\n\n");
			}
		}
	}
	
	//trying to turn the bytes array back into a double
	//double testTwo = *(double*)bytes;
	//printf("byte array back into double = %.3lf\n", testTwo);
	//printf("test + testTwo = %.3lf\n", testTwo+test);
	
	//print results
	printf("Note: single errors turn into double errors, so \nsingle error injection count = single errors injected - double errors injected\n\n");
	printf("Out of %d iterations:\n", iters);
	printf("Single Errors injected: %d\n", se_faults-ded_faults);
	printf("Single Error corrections: %d\n", se_corrections);
	printf("Double Error detections: %d\n", ded_count);
	printf("Double Errors injected: %d\n", ded_faults);
	printf("No Errors found: %d\n", no_error_count);
	
	return 0;
}
*/
