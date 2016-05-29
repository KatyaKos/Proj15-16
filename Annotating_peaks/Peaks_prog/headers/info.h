#ifndef _INFO_H_
#define _INFO_H_

#include  <bits/stdc++.h>
#include "antibody.h"

using namespace std;

extern ifstream fin;
extern ofstream fout;

extern unordered_map<char, long double> am_wght;

class Reader{
private:
	void Read_weights();
	void Read_peaks();
public:
	void Read();
};

class Printer{
private:
	Antibody ant;
	vector<int> mod_pos;

	string lost_atoms(int h, int n, int c);
	void print_seg_peak(const Atom& at, int n);
	void print_mod_peak(const Atom& atl, const Atom& ath);
	void print_pict_peak(const Atom& at, const string& chain);

	void connect(bool (*comp)(const ModifiedChains&, const ModifiedChains&));

	void where_is_cyst(int ci);
	void lonely_cyst(int ci);
	pair<int, int> LRtest(int ci);

	void print_cyst_group(int nC);

public:
	Printer(Antibody ant): ant(ant) {
		mod_pos.assign(peaks.size(), -1);
	}

	void Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&));
	void Pict_Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&));
	void Segments_Cover(bool (*comp)(const ModifiedChains&, const ModifiedChains&));
	void Modified_Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&));
};

#endif