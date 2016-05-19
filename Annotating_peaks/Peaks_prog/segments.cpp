#include "headers/headers.h"

using namespace std;

int peak_search(um_lda& annot, um_ldld& anmass, vector<int>& done, long double delta, const string& chain, int pos){

	int j = 0, num = 0;

	int n = chain.size();
	long double mass = delta, prevmass = delta; delta = 0.0;
	int cyst = 0;

	for(int i = pos; i < n; i++){

		while (done[j]) j++;

		if (i == n - 1) delta = MH2O;

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
			if (check(mass + delta, j)){
				annot[peak] = Atom(pos + 1, i + 1, 0, 0, ic);
				anmass[peak] = mass;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass + delta - MH2O, j)){
				annot[peak] = Atom(pos + 1, i + 1, 1, 0, ic);
				anmass[peak] = mass - MH2O;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass + delta - 2 * MH2O, j)){
				annot[peak] = Atom(pos + 1, i + 1, 2, 0, ic);
				anmass[peak] = mass  - 2 * MH2O;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass + delta - MNH3, j)){
				annot[peak] = Atom(pos + 1, i + 1, 0, 1, ic);
				anmass[peak] = mass - MNH3;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass + delta - MH2O - MNH3, j)){
				annot[peak] = Atom(pos + 1, i + 1, 1, 1, ic);
				anmass[peak] = mass - MH2O - MNH3;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}
			mass -= MCYST * ic;
		}

		if (!flag && mass  > peaks[j] * 1.00001){
			j++;
			mass = prevmass;
			i--;
			if (ch == 'C') cyst--;
		}
	}

	return num;
}