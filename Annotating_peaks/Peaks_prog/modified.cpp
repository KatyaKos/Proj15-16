#include "headers/headers.h"

using namespace std;

int C_preprocess(const string& chain, vector<int>& place){
	int num = 0;
	forn(i, chain.size())
		if (chain[i] == 'C'){
			place.push_back(i);
			num++;
		}
	return num;
}

int check(long double m, int j){
	long double peak = peaks[j];

	if (m <= peak * 1.00001 && m >= peak * 0.99999) return 1;
	return 0;
}

int light_peak_search(long double delta, const string& lchain, int pos, int j){

	int n = lchain.size();
	long double mass = delta + MH2O, prevmass = delta + MH2O;
	int cyst = -1;
	long double peak = peaks[j];

	while (mass + am_wght[lchain[pos]] < peak * 0.99999){
		if (lchain[pos] == 'C'){
			cyst++;
		}
		mass += am_wght[lchain[pos]];
		pos++;

		if (pos == lchain.size()) return 0;
	}

	for(int i = pos; i < n; i++){

		char ch = lchain[i];
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		forn(ic, cyst / 2 + 1){

			mass += MCYST * ic;
			if (check(mass, j)){
				mod_light[peak] = Atom(i + 1, n, 0, 0, ic);
				mass_mod[peak] = mass;
				return 1;
			}else if (check(mass - MH2O, j)){
				mod_light[peak] = Atom(i + 1, n, 1, 0, ic);
				mass_mod[peak] = mass - MH2O;
				break;
			}else if (check(mass - 2 * MH2O, j)){
				mod_light[peak] = Atom(i + 1, n, 2, 0, ic);
				mass_mod[peak] = mass - 2 * MH2O;
				return 1;
			}else if (check(mass - MNH3, j)){
				mod_light[peak] = Atom(i + 1, n, 0, 1, ic);
				mass_mod[peak] = mass - MNH3;
				return 1;
			}else if (check(mass - MH2O - MNH3, j)){
				mod_light[peak] = Atom(i + 1, n, 1, 1, ic);
				mass_mod[peak] = mass - MH2O - MNH3;
				return 1;
			}
			mass -= MCYST * ic;
		}

		if (mass > peak * 1.00001){
			return 0;
		}
	}

	return 0;
}

int modified_peak_search(vector<int>& done, long double delta, int cyst, const string& lchain, const string& hchain, int pos, int posC){

	int j = 0, num = 0;

	int n = hchain.size();
	long double mass = delta, prevmass; delta = 0.0;

	for(int i = posC; i < n; i++){

		while (done[j]) j++;

		if (i == n - 1) delta = MH2O;

		char ch = hchain[i];
		prevmass = mass;
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		int flag = 0;
		
		forn(ic, cyst / 2 + 1){
			long double peak = peaks[j];

			mass += MCYST * ic;
			if (light_peak_search(mass + delta, lchain, 0, j)){
				mod_heavy[peak] = Atom(pos + 1, i + 1, 0, 0, ic);
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MH2O, lchain, 0, j)){
				mod_heavy[peak] = Atom(pos + 1, i + 1, 1, 0, ic);
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (light_peak_search(mass + delta - 2 * MH2O, lchain, 0, j)){
				mod_heavy[peak] = Atom(pos + 1, i + 1, 2, 0, ic);
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MNH3, lchain, 0, j)){
				mod_heavy[peak] = Atom(pos + 1, i + 1, 0, 1, ic);
				done[j] = 1; flag = 1; num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MH2O - MNH3, lchain, 0, j)){
				mod_heavy[peak] = Atom(pos + 1, i + 1, 1, 1, ic);
				done[j] = 1; flag = 1; num++; j++;
				break;
			}
			mass -= MCYST * ic;
		}

		if (!flag && mass > peaks[j] * 1.00001){
			j++;
			mass = prevmass;
			i--;
			if (ch == 'C') cyst--;
		}
	}

	return num;
}