#include "../headers/headers.h"

using namespace std;

int Chains::check(long double m, int j){
	long double peak = peaks[j];

	if (m <= peak * (1 + PREC) && m >= peak * (1 - PREC)) return 1;
	return 0;
}

void Chains::peak_search(um_lda& seg, vector<int>& done, const string& chain, long double mass, int pos){

	int j = 0, nsz = peaks.size();

	int n = chain.size();
	long double prevmass = mass;
	int cyst = 0;

	for(int i = pos; i < n; i++){

		while (done[j]){
			j++;
			if (j == nsz)
				return;
		}

		if (i == n - 1) mass += MH2O;

		char ch = chain[i];
		prevmass = mass;
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		int flag = 0;
		
		forn(ic, cyst / 2 + 1){
			long double peak = peaks[j];

			mass += MCYST * ic;
			if (check(mass, j)){
				seg[peak] = Atom(pos + 1, i + 1, 0, 0, ic);
				mass_seg[peak] = mass;
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (check(mass - MH2O, j)){
				seg[peak] = Atom(pos + 1, i + 1, 1, 0, ic);
				mass_seg[peak] = mass - MH2O;
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (check(mass - 2 * MH2O, j)){
				seg[peak] = Atom(pos + 1, i + 1, 2, 0, ic);
				mass_seg[peak] = mass  - 2 * MH2O;
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (check(mass - MNH3, j)){
				seg[peak] = Atom(pos + 1, i + 1, 0, 1, ic);
				mass_seg[peak] = mass - MNH3;
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (check(mass - MH2O - MNH3, j)){
				seg[peak] = Atom(pos + 1, i + 1, 1, 1, ic);
				mass_seg[peak] = mass - MH2O - MNH3;
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}
			mass -= MCYST * ic;
		}

		if (!flag && mass  > peaks[j] * (1 + PREC)){
			j++;
			mass = prevmass;
			i--;
			if (ch == 'C') cyst--;
		}
	}
}