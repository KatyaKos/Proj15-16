#include "headers/headers.h"

using namespace std;

void LonelyCyst(int ci){

	int i = posC_heavy[ci], psz = posC_heavy.size();
	Printer pr;
	int n = 0;

	for (um_lda::iterator it = mod_heavy.begin(); it != mod_heavy.end(); it++){
		Atom ath = (*it).second, atl = mod_light[(*it).first];
		int flag = 0, lpos = ath.seg.first, rpos = ath.seg.second;

		if (lpos <= i && rpos >= i){
			for (int j = 0; j < psz; j++)
				if (posC_heavy[j] <= rpos && posC_heavy[j] >= lpos && j != ci){
					flag = 1;
					break;
				}

			if (flag) continue;
			string lost = pr.lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ath.numCyst + atl.numCyst);

			fout << setw(36) << right << "ly" + to_string(atl.seg.first) + "+S-S+hi(" + to_string(ath.seg.first) + '-' + to_string(ath.seg.second) + ')' + lost << endl;
			n++;
		}
	}
	fout << endl << "IN TOTAL: " << n << endl;

}

pair<int, int> LRtest(int ci){
	int i = posC_heavy[ci], psz = posC_heavy.size();
	int nseg = 0, n = 0, no = 0;

	for (um_lda::iterator it = mod_heavy.begin(); it != mod_heavy.end(); it++){
		int flag1 = 1, flag2 = 1;
		Atom ath = (*it).second, atl = mod_light[(*it).first];
		int lpos = ath.seg.first, rpos = ath.seg.second;

		if (lpos <= i && rpos >= i){
			n++;
			for (int j = 0; j < psz; j++){
				int jpos = posC_heavy[j];
				if (jpos < i && jpos >= lpos) flag1 = 0;	
				if (j != ci && jpos >= lpos && jpos <= rpos) flag2 = 0;

			}
			if (flag1) nseg++;
			if (flag2) no++;
		}
	}

	if (ci == 0 || ci == psz - 1) return make_pair(n, no);
	return make_pair(nseg, -1);
}