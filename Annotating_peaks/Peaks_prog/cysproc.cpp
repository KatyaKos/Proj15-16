#include "headers.h"

using namespace std;

void fifthCyst(){

	int i = posC_heavy[4];
	Printer pr;
	int n = 0;

	for (um_lda::iterator it = mod_heavy.begin(); it != mod_heavy.end(); it++){
		Atom ath = (*it).second, atl = mod_light[(*it).first];

		if (ath.numCyst == 1 && ath.seg.first <= i && ath.seg.second >= i){
			string lost = pr.lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ath.numCyst + atl.numCyst);

			fout << setw(36) << right << "ly" + to_string(atl.seg.first) + "+S-S+hi(" + to_string(ath.seg.first) + '-' + to_string(ath.seg.second) + ')' + lost << endl;
			n++;
		}
	}
	fout << endl << "IN TOTAL: " << n << endl;

}

int LRtest(int i){
	int nseg = 0;

	for (um_lda::iterator it = mod_heavy.begin(); it != mod_heavy.end(); it++){
		int flag = 0;
		Atom ath = (*it).second, atl = mod_light[(*it).first];
		int lpos = ath.seg.first, rpos = ath.seg.second;

		if (ath.numCyst >= 1 && lpos <= i && rpos >= i){
			for (int j = 0; j < posC_heavy.size(); j++){

				int jpos = posC_heavy[j];
				if (jpos != i && jpos < i && jpos >= lpos){
					flag = 1; break;
				} 
			}
			if (!flag) nseg++;
		}
	}

	return nseg;
}

pair<int, int> FlastSeg(int i){
	int n = 0, no = 0;

	for (um_lda::iterator it = mod_heavy.begin(); it != mod_heavy.end(); it++){
		int flag = 0;
		Atom ath = (*it).second, atl = mod_light[(*it).first];
		int lpos = ath.seg.first, rpos = ath.seg.second;

		if (ath.numCyst >= 1 && lpos <= i && rpos >= i){
			n++;
			if (ath.numCyst == 1) no++;
		}
	}

	return make_pair(n, no);
}