// BranchPredictor.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;

#define SIZE 1048
#define LOCALSIZE 128

int main() {
	string nextLine;
	string addressString;
	long long unsigned int address;
	long long unsigned int previousAddress;
	long long unsigned int progress = 0;
	long long unsigned int hits[5] = { 0, 0, 0, 0, 0 };
	long long unsigned int misses[5] = { 0, 0, 0, 0, 0 };
	long long unsigned int ghits[5] = { 0, 0, 0, 0, 0 };
	long long unsigned int gmisses[5] = { 0, 0, 0, 0, 0 };
	long long unsigned int temp[6] = { 0, 0, 0, 0, 0, 0 };
	long long unsigned int l2hits[5] = { 0, 0, 0, 0 ,0 };
	long long unsigned int l2misses[5] = { 0, 0, 0, 0, 0 };
	long long unsigned int globalRegisters[SIZE*4];
	long long unsigned int gLocalRegisters[SIZE];
	int shiftRegister = 0;
	int shiftRegisterKBits;
	int maxSize;
	int currentN;
	int currentRate;
	bool compulsory = false;
	int updateProgress = 1;
	int FIFO = 0;
	int gFIFO = 0;
	int lFIFO = 0;
	int nCounterBitSize;
	int currentAddress;
	int gCurrentAddress;
	int localShiftRegister;
	char result;
	char prediction;
	char g1Prediction;
	char g2Prediction;
	char gsPrediction;
	char l2Prediction;
	long long unsigned int PHT2[6][SIZE];
	long long unsigned int PHT[6][SIZE];
	long long unsigned int gSHT[6][SIZE];
	long long unsigned int GHT[SIZE];
	int index;
	int gindex;
	int localIndex;
	int local2Address;
	int original;


	int bits = 10;
	int nBits = (pow(2, bits) - 1);

	int testPerformed = 1;

	int assocLMRUResult[5][5];
	int assocFIFOResult[5][5];
	int directMapResult[5][5];

	int l2assocLMRUResult[5][5];
	int l2assocFIFOResult[5][5];
	int l2directMapResult[5][5];

	int l2index;
	int regIndex;

	int g1Hits[5] = { 0, 0, 0, 0, 0 };
	int g2Hits[5] = { 0, 0, 0, 0, 0 };
	int g1Misses[5] = { 0, 0, 0, 0, 0 };
	int g2Misses[5] = { 0, 0, 0, 0, 0 };

	bool local2gShare = false;

	bool fullyAssocLMRU = true;
	bool fullyAssocFIFO = false;
	bool directMap = false;
	bool local2 = true;

	bool g1 = true;
	bool g2 = true;
	bool gShare = true;

	cout << "Please be patient, this will take a few minutes." << endl;
	int C = 0;
	while (C != 2) {
		nBits = (pow(2, bits) - 1);
		//		for (int i = 0; i < (10-bits); i++) {
		//			nBits << 1;
		//		}
		nBits = nBits & 0b11111;
		 fullyAssocLMRU = true;
		 fullyAssocFIFO = false;
		 directMap = false;
		 local2 = true;

		 g1 = true;
		 g2 = true;
		 gShare = true;

		while (testPerformed <= 3) {

			shiftRegister = 0;
			localShiftRegister = 0;

			for (int i = 0; i < SIZE; i++)
			{
				GHT[i] = 0;
			}

			for (int i = 0; i < SIZE; i++) {
				for (int j = 0; j < 6; j++) {
					PHT[j][i] = 0;
					gSHT[j][i] = 0;
					PHT2[j][i] = 0;
				}
			}


			for (int i = 0; i < 5; i++) {
				hits[i] = 0;
				misses[i] = 0;
				g1Hits[i] = 0;
				g2Hits[i] = 0;
				g1Misses[i] = 0;
				g2Misses[i] = 0;
				ghits[i] = 0;
				gmisses[i] = 0;
				l2hits[i] = 0;
				l2misses[i] = 0;
			}

			for (int i = 0; i < 4096; i++) {
				globalRegisters[i] = 0;
			}

			ifstream infile;
			infile.open("trace.din");

			while (getline(infile, nextLine)) {
				for (int i = 0; i < nextLine.length(); i++) {
					if (nextLine[i] == ' ') {
						addressString = nextLine.substr(0, i);
						address = stoull(addressString);
						result = nextLine[i + 1];
						local2Address = address;
						regIndex = address % SIZE*4;
						localIndex = address % LOCALSIZE;
						original = address;

					}
				}


				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (directMap) {
					///////////////////////////////////////////////////////////////
					index = address % SIZE;
					if (PHT[0][index] != 0) {
						PHT[1][index] = 0;
						PHT[2][index] = 0;
						PHT[3][index] = 0;
						PHT[4][index] = 0;
						PHT[5][index] = 0;
					}
					currentAddress = index; //


					if (local2) {
						l2index = ((gLocalRegisters[localIndex] ^ address) & nBits) % SIZE;
						if (local2gShare) address = l2index;
						if (PHT2[0][l2index] != 0) {
							PHT2[1][l2index] = 0;
							PHT2[2][l2index] = 0;
							PHT2[3][l2index] = 0;
							PHT2[4][l2index] = 0;
							PHT2[5][l2index] = 0;
						}
						local2Address = l2index; //
					}

					if (gShare) {
						gindex = ((globalRegisters[regIndex] ^ address) & nBits) % SIZE;
						if (gSHT[0][gindex] != 0) {
							gSHT[1][gindex] = 0;
							gSHT[2][gindex] = 0;
							gSHT[3][gindex] = 0;
							gSHT[4][gindex] = 0;
							gSHT[5][gindex] = 0;
						}
						gCurrentAddress = gindex; //
					}
				}
			
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (fullyAssocFIFO) {
					///////////////////////////////////////////////////////////////
					for (int i = 0; i < SIZE; i++) {
						if (PHT[0][i] == 0) {
							PHT[0][i] = address;
							currentAddress = i;
							break;
						}
						//////////////////////////////
						else if (PHT[0][i] == address) {
							currentAddress = i;
							break;
						}
						//////////////////////////////
						else if (i == (SIZE-1)) {
							if (FIFO == SIZE) FIFO = 0;
							PHT[0][FIFO] = address;
							PHT[0][FIFO] = 0;
							currentAddress = FIFO;
							FIFO++;
						}
					}
					if (local2) {
						local2Address = ((original ^ gLocalRegisters[localIndex]) & nBits);
						if (local2gShare) address = local2Address;
						for (int i = 0; i < SIZE; i++) {
							if (PHT2[0][i] == 0) {
								PHT2[0][i] = address;
								local2Address = i;
								break;
							}
							//////////////////////////////
							else if (gSHT[0][i] == local2Address) {
								local2Address = i;
								break;
							}
							//////////////////////////////
							else if (i == (SIZE-1)) {
								if (lFIFO == SIZE) lFIFO = 0;
								gSHT[0][lFIFO] = local2Address;
								gSHT[0][lFIFO] = 0;
								local2Address = lFIFO;
								gFIFO++;
							}
						}
					}
					if (gShare) {
						address = ((address ^ globalRegisters[regIndex]) & nBits);
						for (int i = 0; i < SIZE; i++) {
							if (gSHT[0][i] == 0) {
								gSHT[0][i] = address;
								gCurrentAddress = i;
								break;
							}
							//////////////////////////////
							else if (gSHT[0][i] == address) {
								gCurrentAddress = i;
								break;
							}
							//////////////////////////////
							else if (i == (SIZE-1)) {
								if (gFIFO == SIZE) gFIFO = 0;
								gSHT[0][gFIFO] = address;
								gSHT[0][gFIFO] = 0;
								gCurrentAddress = gFIFO;
								gFIFO++;
							}
						}
					}
				}
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (fullyAssocLMRU) {
					///////////////////////////////////////////////////////////////
					for (int i = 0; i < SIZE; i++) {
						if (PHT[0][i] == 0) {
							for (int j = 0; j < i; j++) {
								for (int k = 0; k < 6; k++) {
									PHT[k][i - j] = PHT[k][i - 1 - j];
								}
							}
							for (int j = 0; j < 6; j++) {
								PHT[j][0] = 0;
							}
							PHT[0][0] = address;
							currentAddress = 0;
							break;
						}
						//////////////////////////////
						else if (PHT[0][i] == address) {
							for (int j = 0; j < 6; j++) {
								temp[j] = PHT[j][i];
							}
							//////////////////////////////
							for (int j = 0; j < i; j++) {
								for (int k = 0; k < 6; k++) {
									PHT[k][i - j] = PHT[k][i - 1 - j];
								}
							}
							//////////////////////////////
							for (int j = 0; j < 6; j++) {
								PHT[j][0] = temp[j];
								currentAddress = 0;
							}
							break;
						}
						//////////////////////////////
						else if (i == (SIZE-1)) {
							for (int j = 0; j < (SIZE-1); j++) {
								for (int k = 0; k < 6; k++) {
									PHT[k][(SIZE-1) - j] = PHT[k][(SIZE-2) - j];
								}
							}
							//////////////////////////////
							for (int j = 0; j < 6; j++) {
								PHT[j][0] = 0;
							}
							PHT[0][0] = address;
							currentAddress = 0;
							break;
						}
					}
					if (local2) {
						local2Address = ((original ^ gLocalRegisters[localIndex]) & nBits);
						if (local2gShare) address = local2Address;
						for (int i = 0; i < SIZE; i++) {
							if (PHT2[0][i] == 0) {
								for (int j = 0; j < i; j++) {
									for (int k = 0; k < 6; k++) {
										PHT2[k][i - j] = PHT2[k][i - 1 - j];
									}
								}
								for (int j = 0; j < 6; j++) {
									PHT2[j][0] = 0;
								}
								PHT2[0][0] = local2Address;
								local2Address = 0;
								break;
							}
							//////////////////////////////
							else if (PHT2[0][i] == local2Address) {
								for (int j = 0; j < 6; j++) {
									temp[j] = PHT2[j][i];
								}
								//////////////////////////////
								for (int j = 0; j < i; j++) {
									for (int k = 0; k < 6; k++) {
										PHT2[k][i - j] = PHT2[k][i - 1 - j];
									}
								}
								//////////////////////////////
								for (int j = 0; j < 6; j++) {
									PHT2[j][0] = temp[j];
								}
								local2Address = 0;
								break;
							}
							//////////////////////////////
							else if (i == (SIZE-1)) {
								for (int j = 0; j < (SIZE-1); j++) {
									for (int k = 0; k < 6; k++) {
										PHT2[k][(SIZE-1) - j] = PHT2[k][(SIZE-2) - j];
									}
								}
								//////////////////////////////
								for (int j = 0; j < 6; j++) {
									PHT2[j][0] = 0;
								}
								PHT2[0][0] = local2Address;
								local2Address = 0;
								break;
							}
						}
					}
					if (gShare) {
						address = ((address ^ globalRegisters[regIndex]) & nBits);
						for (int i = 0; i < SIZE; i++) {
							if (gSHT[0][i] == 0) {
								for (int j = 0; j < i; j++) {
									for (int k = 0; k < 6; k++) {
										gSHT[k][i - j] = gSHT[k][i - 1 - j];
									}
								}
								for (int j = 0; j < 6; j++) {
									gSHT[j][0] = 0;
								}
								gSHT[0][0] = address;
								gCurrentAddress = 0;
								break;
							}
							//////////////////////////////
							else if (gSHT[0][i] == address) {
								for (int j = 0; j < 6; j++) {
									temp[j] = gSHT[j][i];
								}
								//////////////////////////////
								for (int j = 0; j < i; j++) {
									for (int k = 0; k < 6; k++) {
										gSHT[k][i - j] = gSHT[k][i - 1 - j];
									}
								}
								//////////////////////////////
								for (int j = 0; j < 6; j++) {
									gSHT[j][0] = temp[j];
								}
								gCurrentAddress = 0;
								break;
							}
							//////////////////////////////
							else if (i == (SIZE-1)) {
								for (int j = 0; j < (SIZE-1); j++) {
									for (int k = 0; k < 6; k++) {
										gSHT[k][(SIZE-1) - j] = gSHT[k][(SIZE-2) - j];
									}
								}
								//////////////////////////////
								for (int j = 0; j < 6; j++) {
									gSHT[j][0] = 0;
								}
								gSHT[0][0] = address;
								gCurrentAddress = 0;
								break;
							}
						}
					}
				}

				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				nCounterBitSize = 2;
				shiftRegister = shiftRegister & 0b1111111111;
				localShiftRegister = localShiftRegister & 0b1111111111;
				shiftRegisterKBits = 0;

				while (nCounterBitSize <= 10) {
					maxSize = (pow(2, nCounterBitSize) - 1);
					shiftRegisterKBits = shiftRegister & maxSize;
					currentN = nCounterBitSize / 2;
					currentRate = currentN - 1;
					///////////////////////////////////////////////////////////////
					if (PHT[currentN][currentAddress] <= (maxSize / (2.0))) prediction = 'N';
					else prediction = 'T';
					///////////////////////////////////////////////////////////////
					if (prediction == result) {
						hits[currentRate]++;
						if (prediction == 'T') {
							if (PHT[currentN][currentAddress] < (maxSize - 1)) PHT[currentN][currentAddress]++;
						}
						if (prediction == 'N') {
							if (PHT[currentN][currentAddress] > 0) PHT[currentN][currentAddress]--;
						}
					}
					//////////////////////////////
					else {
						misses[currentRate]++;
						if (prediction == 'N') {
							if (PHT[currentN][currentAddress] < (maxSize - 1)) PHT[currentN][currentAddress]++;
						}
						if (prediction == 'T') {
							if (PHT[currentN][currentAddress] > 0) PHT[currentN][currentAddress]--;
						}
					}
					///////////////////////////////////////////////////////////
					if (gSHT[currentN][gCurrentAddress] <= (maxSize / (2.0))) gsPrediction = 'N';
					else gsPrediction = 'T';
					///////////////////////////////////////////////////////////////
					if (gsPrediction == result) {
						ghits[currentRate]++;
						if (prediction == 'T') {
							if (gSHT[currentN][gCurrentAddress] < (maxSize - 1)) gSHT[currentN][gCurrentAddress]++;
						}
						if (gsPrediction == 'N') {
							if (gSHT[currentN][gCurrentAddress] > 0) gSHT[currentN][gCurrentAddress]--;
						}
					}
					//////////////////////////////
					else {
						gmisses[currentRate]++;
						if (gsPrediction == 'N') {
							if (gSHT[currentN][gCurrentAddress] < (maxSize - 1)) gSHT[currentN][gCurrentAddress]++;
						}
						if (gsPrediction == 'T') {
							if (gSHT[currentN][gCurrentAddress] > 0) gSHT[currentN][gCurrentAddress]--;
						}
					}
					///////////////////////////////////////////////////////////
					//1-Level Global
					if (g1) {
						if (shiftRegisterKBits <= (maxSize / (2.0))) g1Prediction = 'N';
						else g1Prediction = 'T';
						//////////////////////////////
						if (g1Prediction == result) g1Hits[currentRate]++;
						else g1Misses[currentRate]++;
					}
					///////////////////////////////////////////////////////////
					//2-Level Global
					if (g2) {
						if (GHT[shiftRegister] <= (maxSize / (2.0))) g2Prediction = 'N';
						else g2Prediction = 'T';

						if (g2Prediction == result) {
							g2Hits[currentRate]++;
							if (g2Prediction == 'T') {
								if (GHT[shiftRegister] < (maxSize - 1)) GHT[shiftRegister]++;
							}
							if (g2Prediction == 'N') {
								if (GHT[shiftRegister] > 0) GHT[shiftRegister]--;
							}
						}
						//////////////////////////////
						else {
							g2Misses[currentRate]++;
							if (g2Prediction == 'N') {
								if (GHT[shiftRegister] < (maxSize - 1)) GHT[shiftRegister]++;
							}
							if (g2Prediction == 'T') {
								if (GHT[shiftRegister] > 0) GHT[shiftRegister]--;
							}

						}

					}
					if (local2) {
						if (PHT2[currentN][localShiftRegister] <= (maxSize / (2.0))) l2Prediction = 'N';
						else l2Prediction = 'T';

						if (l2Prediction == result) {
							l2hits[currentRate]++;
							if (l2Prediction == 'T') {
								if (PHT2[currentN][localShiftRegister] < (maxSize - 1)) PHT2[currentN][localShiftRegister]++;
							}
							if (l2Prediction == 'N') {
								if (PHT2[currentN][localShiftRegister] > 0) PHT2[currentN][localShiftRegister]--;
							}
						}
						//////////////////////////////
						else {
							l2misses[currentRate]++;
							if (l2Prediction == 'N') {
								if (PHT2[currentN][localShiftRegister] < (maxSize - 1)) PHT2[currentN][localShiftRegister]++;
							}
							if (l2Prediction == 'T') {
								if (PHT2[currentN][shiftRegister] > 0) PHT2[currentN][shiftRegister]--;
							}

						}

					}

					///////////////////////////////////////////////////////////////
					//gShare
					if (gShare) {
						if (result == 'T') globalRegisters[regIndex] = ((globalRegisters[regIndex] << 1) + 1) & 0b1111111111;
						else globalRegisters[regIndex] = ((globalRegisters[regIndex] << 1)) & 0b1111111111;
					}


					///////////////////////////////////////////////////////////////
					nCounterBitSize += 2;
				} //End of result calculation while loop
				///////////////////////////////////////////////////////////////
				if (result == 'T') (shiftRegister << 1) + 1;
				else shiftRegister << 1;

				if (result == 'T') (localShiftRegister << 1) + 1;
				else localShiftRegister << 1;

				progress++;
				if (progress >= 164162) {
					if (updateProgress == 100) updateProgress = 0;
					cout << "\r" << "Test #" << testPerformed << " Progress: " << updateProgress << "% Complete. Round " << C+1;
					updateProgress++;
					progress = 0;
				}

			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (fullyAssocLMRU) {
				for (int i = 0; i < 5; i++) {
					if (!local2gShare) {
						assocLMRUResult[0][i] = (i + 1) * 2;
						assocLMRUResult[1][i] = hits[i];
						assocLMRUResult[2][i] = misses[i];
						assocLMRUResult[3][i] = ghits[i];
						assocLMRUResult[4][i] = gmisses[i];
					}

					l2assocLMRUResult[0][i] = (i + 1) * 2;
					if (!local2gShare) l2assocLMRUResult[1][i] = l2hits[i];
					if (!local2gShare) l2assocLMRUResult[2][i] = l2misses[i];
					if (local2gShare) l2assocLMRUResult[3][i] = ghits[i];
					if (local2gShare) l2assocLMRUResult[4][i] = gmisses[i];
				}
				fullyAssocLMRU = false;
				fullyAssocFIFO = true;
			}
			///////////////////////////////////////////////////////////////
			else if (fullyAssocFIFO) {
				for (int i = 0; i < 5; i++) {
					if (!local2gShare) {
						assocFIFOResult[0][i] = (i + 1) * 2;
						assocFIFOResult[1][i] = hits[i];
						assocFIFOResult[2][i] = misses[i];
						assocFIFOResult[3][i] = ghits[i];
						assocFIFOResult[4][i] = gmisses[i];
					}

					l2assocFIFOResult[0][i] = (i + 1) * 2;
					if (!local2gShare) l2assocFIFOResult[1][i] = l2hits[i];
					if (!local2gShare) l2assocFIFOResult[2][i] = l2misses[i];
					if (local2gShare) l2assocFIFOResult[3][i] = ghits[i];
					if (local2gShare) l2assocFIFOResult[4][i] = gmisses[i];
				}
				fullyAssocFIFO = false;
				directMap = true;
			}
			///////////////////////////////////////////////////////////////
			else if (directMap) {
				for (int i = 0; i < 5; i++) {
					if (!local2gShare) {
						directMapResult[0][i] = (i + 1) * 2;
						directMapResult[1][i] = hits[i];
						directMapResult[2][i] = misses[i];
						directMapResult[3][i] = ghits[i];
						directMapResult[4][i] = gmisses[i];
					}

					l2directMapResult[0][i] = (i + 1) * 2;
					if (!local2gShare) l2directMapResult[1][i] = l2hits[i];
					if (!local2gShare) l2directMapResult[2][i] = l2misses[i];
					if (local2gShare) l2directMapResult[3][i] = ghits[i];
					if (local2gShare) l2directMapResult[4][i] = gmisses[i];
				}
				directMap = false;
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			infile.close();
			testPerformed++;
		}

		testPerformed = 1;
		local2gShare = true;
		C++;
	}

		ofstream cout("output.txt");

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cout << endl;
		cout << "Branch Prediction Results: " << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "Fully Associative: LMRU Replacement Policy" << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		for (int i = 0; i < 5; i++) {
			cout << endl << "Local:                " << assocLMRUResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((assocLMRUResult[1][i] * 1.0) / (assocLMRUResult[1][i] + assocLMRUResult[2][i]) * 100) << endl;
			cout << endl << "gShare:               " << assocLMRUResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((assocLMRUResult[3][i] * 1.0) / (assocLMRUResult[3][i] + assocLMRUResult[4][i]) * 100) << endl;
			cout << endl << "Level 2 Local:        " << l2assocLMRUResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2assocLMRUResult[1][i] * 1.0) / (l2assocLMRUResult[1][i] + l2assocLMRUResult[2][i]) * 100) << endl;
			cout << endl << "Level 2 Local gShare: " << l2assocLMRUResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2assocLMRUResult[3][i] * 1.0) / (l2assocLMRUResult[3][i] + l2assocLMRUResult[4][i]) * 100) << endl;
			cout << "---------------------------------------------------------------" << endl;
		}
		cout << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "Fully Associative: FIFO Replacement Policy" << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		for (int i = 0; i < 5; i++) {
			cout << endl << "Local:                " << assocFIFOResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((assocFIFOResult[1][i] * 1.0) / (assocFIFOResult[1][i] + assocFIFOResult[2][i]) * 100) << endl;
			cout << endl << "gShare:               " << assocFIFOResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((assocFIFOResult[3][i] * 1.0) / (assocFIFOResult[3][i] + assocFIFOResult[4][i]) * 100) << endl;
			cout << endl << "Level 2 Local:        " << l2assocFIFOResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2assocFIFOResult[1][i] * 1.0) / (l2assocFIFOResult[1][i] + l2assocFIFOResult[2][i]) * 100) << endl;
			cout << endl << "Level 2 Local gShare: " << l2assocFIFOResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2assocFIFOResult[3][i] * 1.0) / (l2assocFIFOResult[3][i] + l2assocFIFOResult[4][i]) * 100) << endl;
			cout << "---------------------------------------------------------------" << endl;
		}
		cout << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "Direct Mapped:" << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		for (int i = 0; i < 5; i++) {
			cout << endl << "Local:                " << directMapResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((directMapResult[1][i] * 1.0) / (directMapResult[1][i] + directMapResult[2][i]) * 100) << endl;
			cout << endl << "gShare:               " << directMapResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((directMapResult[3][i] * 1.0) / (directMapResult[3][i] + directMapResult[4][i]) * 100) << endl;
			cout << endl << "Level 2 Local:        " << l2directMapResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2directMapResult[1][i] * 1.0) / (l2directMapResult[1][i] + l2directMapResult[2][i]) * 100) << endl;
			cout << endl << "Level 2 Local gShare: " << l2directMapResult[0][i] << "-bit Counter: " << "Prediction Rate = " << ((l2directMapResult[3][i] * 1.0) / (l2directMapResult[3][i] + l2directMapResult[4][i]) * 100) << endl;
			cout << "---------------------------------------------------------------" << endl;
		}
		cout << endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "One Layer Global Branch Prediction Results: " << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		for (int i = 0; i < 5; i++) {
			cout << endl << (i + 1) * 2 << "-bit Counter: Prediction Rate = " << ((g1Hits[i] * 1.0) / (g1Hits[i] + g1Misses[i])) * 100 << endl;
			cout << "---------------------------------------------------------------" << endl;
		}
		cout << endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cout << "///////////////////////////////////////////////////////////////" << endl;
		cout << "Two Layer Global Branch Prediction Results: " << endl;
		cout << "///////////////////////////////////////////////////////////////" << endl;
		for (int i = 0; i < 5; i++) {
			cout << endl << (i + 1) * 2 << "-bit Counter: Prediction Rate = " << ((g2Hits[i] * 1.0) / (g2Hits[i] + g2Misses[i])) * 100 << endl;
			cout << "---------------------------------------------------------------" << endl;
		}
		return 0;
	}


