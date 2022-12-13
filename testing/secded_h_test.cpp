//james bickerstaff
//testing basic double data type functions for SEC-DED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <cstdlib>
#include "secded.h"

int main(int argc, char **argv)
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
	printf("Please input a number to test: ");
	std::cin >> test;
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
	
	//verifying taking data from char array of 9 back into a double
	double backToDouble = charToDouble(transmit);
	if(backToDouble != test)
		printf("DATA TURNED BACK INTO DOUBLE DOES NOT MATCH INPUT DATA!\n\n");
	else
		printf("DATA BACK TO DOUBLE SUCCESSFULLY!\n\n");
	
	//done testing, output results
	//print results
	printf("Out of %d iterations:\n", iters);
	printf("Single Errors injected: %d\n", se_faults);
	printf("Single Error corrections: %d\n", se_corrections);
	printf("Double Error detections: %d\n", ded_count);
	printf("Double Errors injected: %d\n", ded_faults);
	printf("No Errors found: %d\n", no_error_count);
}