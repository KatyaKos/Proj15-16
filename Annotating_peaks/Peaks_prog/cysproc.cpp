#include "headers.h"

using namespace std;

void fifthCyst(){
	fout.open(CYS_PROCESS_FILE);

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

	fout.close();

}