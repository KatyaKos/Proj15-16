#ifndef _INFO_H_
#define _INFO_H_

#include  <bits/stdc++.h>
#include "antibody.h"

using namespace std;

extern ifstream fin;
extern ofstream fout;

class Spectrum{
private:
	string spect, anti;
	string annot, modif, pict, segm;

	void Read_peaks();

	Antibody ant;
	vector<int> mod_pos;

	string lost_atoms(int h, int n, int c);
	void print_seg_peak(const Atom& at, int n);
	void print_mod_peak(const Atom& atl, const Atom& ath);
	void print_pict_peak(const Atom& at, const string& chain);

	void connect(bool (*comp)(const ModChains&, const ModChains&));

	void where_is_cyst(int ci);
	void lonely_cyst(int ci);
	pair<int, int> LRtest(int ci);

	void print_cyst_group(int nC);
public:
	Spectrum(string s, string a): spect(s), anti(a) {}

	void Read();
	void set_antibody(Antibody anti){
		ant = anti;
		mod_pos.assign(peaks.size(), -1);
	}
	void Annotate(bool (*comp)(const ModChains&, const ModChains&));
	void Pict_Annotate(bool (*comp)(const ModChains&, const ModChains&));
	void Segments_Cover(bool (*comp)(const ModChains&, const ModChains&));
	void Modified_Annotate(bool (*comp)(const ModChains&, const ModChains&));
};

#endif