#ifndef _INFO_H_
#define _INFO_H_

#include  <bits/stdc++.h>

using namespace std;

#define forn(i, n) for (int i = 0; i < n; i++)

struct Atom{
	int numH2O, numNH3, numCyst;
	pair <int, int> seg;\
	Atom() = default;
	Atom(int i, int j, int h, int n, int c): seg(make_pair(i, j)), numCyst(c), numNH3(n), numH2O(h) {}
};

typedef unordered_map<long double, Atom> um_lda;
typedef unordered_map<long double, long double> um_ldld;

extern ifstream fin;
extern ofstream fout;

extern unordered_map<char, long double> am_wght;

extern um_lda seg_light;
extern um_lda seg_heavy;
extern um_lda mod_light;
extern um_lda mod_heavy;

extern um_ldld mass_seg_light;
extern um_ldld mass_seg_heavy;
extern um_ldld mass_mod;

extern vector<int> posC_heavy;
extern vector<long double> peaks;

extern int seg_num_heavy, seg_num_light, mod_num;

#endif