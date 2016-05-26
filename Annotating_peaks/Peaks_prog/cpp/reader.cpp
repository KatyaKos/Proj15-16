#include "../headers/headers.h"

using namespace std;

void Reader::Read(){
	Read_weights();
	Read_peaks();
}

void Reader::Read_weights(){
	fin.open(AMINO_FILE);
	char name, tmp;
	long double w;
	while (fin >> name){
		fin >> tmp >> w >> tmp;
		am_wght[name] = w;
	}
	fin.close();
}

void Reader::Read_peaks(){
	fin.open(SPECTRUM_FILE);
	string buff;
	forn(i, 6){
		getline(fin, buff);
	}
	long double d;

	while (fin >> d){
		peaks.push_back(d);
		fin >> d >> d;
	}

	sort(peaks.begin(), peaks.end());
	fin.close();
}