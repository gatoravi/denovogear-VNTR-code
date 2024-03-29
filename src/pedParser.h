#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#define MAX_TRIOS 1000 // max number of trios in the bcf filed
#define MAX_PAIRS 100 // max number of paired samples in the bcf file
#define ID_LENGTH 1000
#define LINE_LENGTH 10000

// Trio Structure
typedef struct  {
  char fID[ID_LENGTH];  // family ID
  char cID[ID_LENGTH];  // child sampleID
  char dID[ID_LENGTH];  // dad sampleID
  char mID[ID_LENGTH];  // mom sampleID
} Trio;

// Pair Structure
typedef struct  {
  char pairID[ID_LENGTH]; // case ID
  char tumorID[ID_LENGTH];  // tumor sample ID
  char normalID[ID_LENGTH]; // normal sample ID
} Pair;

void parse_ped(std::string ped_file, Trio** t, Pair** p, int& trio_count, int& pair_count);