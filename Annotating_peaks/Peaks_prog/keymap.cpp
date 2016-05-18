#include "headers/headers.h"

using namespace std;

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

void SegCover(){
	fout.open(CYS_PROCESS_FILE);
	fout << "SEGMENTS THAT CONTAIN ONLY FIFTH CYSTEIN:" << endl << endl;
	LonelyCyst(4);

	int psz = posC_heavy.size();
	fout << endl << endl << endl << "NUMBER OF SEGMENTS THAT COVER CYSTEINS WITH THEIR LEFTEST/RIGHTEST END:" << endl;
	for (int i = 1; i < psz - 1; i++){
		fout << "cystein " << i + 1 << ":  " <<  LRtest(i).first << endl;
	}

	pair<int, int> pf = LRtest(0), pl = LRtest(psz - 1);
	fout << endl << endl << "NUMBER OF SEGMENTS THAT CONTAIN Nth CYSTEIN:" << endl;
	fout << "cystein 1:  " << pf.first << endl << "cystein " << psz << ":  " << pl.first << endl << endl;
	fout << endl << endl << "NUMBER OF SEGMENTS THAT CONTAIN ONLY Nth CYSTEIN:" << endl;
	fout << "cystein 1:  " << pf.second << endl << "cystein " << psz << ":  " << pl.second << endl << endl;

	fout.close();
}