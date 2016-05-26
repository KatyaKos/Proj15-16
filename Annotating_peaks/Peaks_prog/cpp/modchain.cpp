#include "../headers/headers.h"

using namespace std;

int ModifiedChains::light_peak_search(long double delta, int pos, int j){

	int n = ant.lchain.size();
	long double mass = delta + MH2O, prevmass = delta + MH2O;
	int cyst = -1;
	long double peak = peaks[j];

	while (mass + am_wght[ant.lchain[pos]] < peak * 0.99999){
		if (ant.lchain[pos] == 'C'){
			cyst++;
		}
		mass += am_wght[ant.lchain[pos]];
		pos++;

		if (pos == ant.lchain.size()) return 0;
	}

	for(int i = pos; i < n; i++){

		char ch = ant.lchain[i];
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		forn(ic, cyst / 2 + 1){

			mass += MCYST * ic;
			if (check(mass, j)){
				seg_light[peak] = Atom(n - i, n, 0, 0, ic);
				mass_seg[peak] = mass;
				return 1;
			}else if (check(mass - MH2O, j)){
				seg_light[peak] = Atom(n - i, n, 1, 0, ic);
				mass_seg[peak] = mass - MH2O;
				return 1;
			}else if (check(mass - 2 * MH2O, j)){
				seg_light[peak] = Atom(n - i, n, 2, 0, ic);
				mass_seg[peak] = mass - 2 * MH2O;
				return 1;
			}else if (check(mass - MNH3, j)){
				seg_light[peak] = Atom(n - i, n, 0, 1, ic);
				mass_seg[peak] = mass - MNH3;
				return 1;
			}else if (check(mass - MH2O - MNH3, j)){
				seg_light[peak] = Atom(n - i, n, 1, 1, ic);
				mass_seg[peak] = mass - MH2O - MNH3;
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

void ModifiedChains::peak_search(long double delta, int cyst, int pos, int nC){

	int j = 0, posC = ant.posC_heavy[nC], nsz = peaks.size();

	int n = ant.hchain.size();
	long double mass = delta, prevmass; delta = 0.0;

	for(int i = posC; i < n; i++){

		while (done[j]){
			j++;
			if (j == nsz) return;
		}

		if (i == n - 1) delta = MH2O;

		char ch = ant.hchain[i];
		prevmass = mass;
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		int flag = 0;
		long double peak = peaks[j];
		
		forn(ic, cyst / 2 + 1){

			mass += MCYST * ic;
			if (light_peak_search(mass + delta, 0, j)){
				seg_heavy[peak] = Atom(pos + 1, i + 1, 0, 0, ic);
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MH2O, 0, j)){
				seg_heavy[peak] = Atom(pos + 1, i + 1, 1, 0, ic);
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (light_peak_search(mass + delta - 2 * MH2O, 0, j)){
				seg_heavy[peak] = Atom(pos + 1, i + 1, 2, 0, ic);
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MNH3, 0, j)){
				seg_heavy[peak] = Atom(pos + 1, i + 1, 0, 1, ic);
				done[j] = 1; flag = 1; seg_num++; j++;
				break;
			}else if (light_peak_search(mass + delta - MH2O - MNH3, 0, j)){
				seg_heavy[peak] = Atom(pos + 1, i + 1, 1, 1, ic);
				done[j] = 1; flag = 1; seg_num++; j++;
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
}

void ModifiedChains::chain_process(int ii){

	int i = ant.posC_heavy[ii], hcyst = -1;
	long double hMass = -am_wght['C'] + MCYST;
	cout << i << endl;
	done.assign(peaks.size(), 0);
	seg_num = 0;

	for (int rpos = i; rpos >= 0; rpos--){
		char ch = ant.hchain[rpos];
		if (ch == 'C') hcyst++;
		hMass += am_wght[ch];
		peak_search(hMass, hcyst - 1, rpos, ii);
	}

}