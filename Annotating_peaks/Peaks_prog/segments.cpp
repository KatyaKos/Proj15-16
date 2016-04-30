#include "headers.h"

using namespace std;

int peak_search(um_lda& annot, um_ldld& anmass, vector<int>& done, long double delta, string chain, int pos){

	int j = 0, num = 0;

	int n = chain.size();
	long double mass = delta, prevmass = delta;
	int cyst = 0;

	for(int i = pos; i < n; i++){

		if (pos && i == n - 1) break;

		while (done[j]) j++;

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
				annot[peak] = Atom(pos + 1, i + 1, 0, 0, ic);
				anmass[peak] = mass;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass - MH2O, j)){
				annot[peak] = Atom(pos + 1, i + 1, 1, 0, ic);
				anmass[peak] = mass - MH2O;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass - 2 * MH2O, j)){
				annot[peak] = Atom(pos + 1, i + 1, 2, 0, ic);
				anmass[peak] = mass  - 2 * MH2O;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass - MNH3, j)){
				annot[peak] = Atom(pos + 1, i + 1, 0, 1, ic);
				anmass[peak] = mass - MNH3;
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (check(mass - MH2O - MNH3, j)){
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