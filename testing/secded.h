#ifndef SECDED_H
#define SECDED_H

//FLOW FOR USAGE OF SEC/DED CODE:
// 1) must seed srand for more random spread of bit flips: 
//    1.1) "std::srand(time(0));"
// 2) also want to take as input from cfg or cmd or something the values for chance of single-bit and double-bit errors
//    2.1) "singleErrorChance" and "doubleErrorChance" in this file
// 3) get double from program and create a 9 byte UNSIGNED char array for it (8 bytes data and 1 byte check)
//    3.1) do not modify the char array, just create: "unsigned char transmit[9];"
// 4) pass the double and 9-char array to the method "addCheckByte(double data, unsigned char codeWord[])"
//    4.1) this will turn the data into a byte array and tack on the check byte as well, just edits the codeWord input directly, no returns
// 5) modify bit(s) within the byte array using "modifyBits(unsigned char codeWord[], int singleErrorChance, int doubleErrorChance)" or other function
//    5.1) this function returns the number of bits modified (0, 1, or 2)
// 6) transmit and receive data with MPI calls
// 7) use the function "checkReceivedData(unsigned char receivedData[])" by passing in the received unsigned char array
//    7.1) this will correct data if possible, or flag a double bit error (any fixes made are applied directly to the receivedData array)
//    7.2) int returned denoting if: no error detected (return=0), 1 error fixed (return=1), or double error detected (return=2)
// 8) turn data back into a double by using the "double charToDouble(unsigned char receivedData[])" method
//    8.1) just pass in the corrected data from step 7, and it will return the data as a double
// 9) might want to count what the error was and if it was fixed to display test data at the end, possible through return values of "modifyBits" and "checkReceivedData"


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

//BEGIN SUPPORTING METHODS - JUMP DOWN TO "addCheckByte" FOR METHODS TO CALL WITHIN PROGRAM

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

//END SUPPORTING METHODS


//BEGIN METHODS TO CALL FROM PROGRAM 

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
// int denoting if: no error detected (0), 1 error fixed (1), or double error detected (2)
// also directly fix the data array of received data
int checkReceivedData(unsigned char receivedData[])
{	
	//bool for if data has been fixed (or not needed fixed) or needs retransmitted (DED)
	int correctedCount = -1;
	//calculate syndrome based on received data
	unsigned char syndrome = calculateSyndrome(receivedData);
	
	//do all the following if syndrome != 0, meaning no error
	if(syndrome != 0x00)
	{
		//correct bit if possible
		bool fixedInfo = secded_data(receivedData, syndrome);
		
		//if it could not be corrected,double error has been detected
		if(fixedInfo)
		{
			//DATA SUCCESSFULLY FIXED :D
			//set the return value to 1, as the data has been fixed
			correctedCount = 1;
		}
		else
		{
			//DED DETECTED, NEED TO RETRANSMIT DATA
			//set return value to 2, as DED has been encountered
			correctedCount = 2;
		}		
	}
	else
	{
		//NO ERROR DETECTED, DATA NEEDS NO MODIFICATION
		//make return value 0, as no errors encountered
		correctedCount = 0;
	}
	
	//return bool signifying if data is good or not
	return correctedCount;
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

//END METHODS TO CALL FROM PROGRAM

#endif