#include "headers/headers.h"

using namespace std;

ifstream fin;
ofstream fout;

unordered_map<char, long double> am_wght;

um_lda seg_light;
um_lda seg_heavy;
um_lda mod_light;
um_lda mod_heavy;

um_ldld mass_seg_light;
um_ldld mass_seg_heavy;
um_ldld mass_mod;

vector<int> posC_heavy;
vector<long double> peaks;

int seg_num_heavy, seg_num_light, mod_num;

int main(){

	Read_weights();
	Read_peaks();

	string lchain, hchain;
	fin.open(CHAINS_FILE);
	getline(fin, lchain);
	getline(fin, lchain);
	getline(fin, hchain);
	getline(fin, hchain);
	fin.close();

	seg_num_light = chains_normal(lchain, seg_light, mass_seg_light);
	seg_num_heavy = chains_normal(hchain, seg_heavy, mass_seg_heavy);
	chains_modified(lchain, hchain);

	annotating(lchain, hchain);
	pict_annotating(lchain, hchain);

	SegCover();

	return 0;
}