#include "headers.h"

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void chains_modified(string lchain, string hchain){

	vector <int> done;
	done.assign(peaks.size(), 0);

	lchain.pop_back();
	reverse(lchain.begin(), lchain.end());
	hchain.pop_back();

	int num_hC = C_preprocess(hchain, posC_heavy);
	mod_num = 0;

	forn(ii, num_hC){

		int i = posC_heavy[ii], hcyst = -1;
		long double hMass = -am_wght['C'] + MCYST;
		cout << i << endl;

		for (int rpos = i; rpos >= 0; rpos--){
			char ch = hchain[rpos];
			if (ch == 'C') hcyst++;
			hMass += am_wght[ch];
			mod_num += modified_peak_search(done, hMass, hcyst - 1, lchain, hchain, rpos, i);
		}

	}

}


int chains_normal(string chain, um_lda& annot, um_ldld& anmass){

	vector <int> done;
	done.assign(peaks.size(), 0);

	chain.pop_back();
	int num = 0;

	for(int i = 1; i < chain.size() - 1; i++){
		num += peak_search(annot, anmass, done, 0.0, chain, i);
	}

	num += peak_search(annot, anmass, done, 0.0, chain, 0);
	reverse(chain.begin(), chain.end());
	num += peak_search(annot, anmass, done, MH2O, chain, 0);

	return num;
}

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

	check_maps(lchain, hchain);

	return 0;
}