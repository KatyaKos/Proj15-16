#include "../headers/headers.h"

using namespace std;

void ModifiedChains::heavy_search(long double mass, int nH, int nN, int nC, int pos){
	string chain = ant.hchain;
	int n = chain.size(), psz = done.size(), nl = ant.lchain.size();

	forn(i, n){
		peak_search(seg_heavy, done, chain, mass, i);
		int l, r;
		forn(j, psz){
			if (done[j] != 1) continue;
			long double peak = peaks[j];
			l = lower_bound(ant.posC_heavy.begin(), ant.posC_heavy.end(), seg_heavy[peak].seg.first - 1) - ant.posC_heavy.begin();
			r = upper_bound(ant.posC_heavy.begin(), ant.posC_heavy.end(), seg_heavy[peak].seg.second - 1) - ant.posC_heavy.begin();
			//cout << l << ' ' << r << endl;

			if (r - l <= 0){
				done[j] = 0;
				mass_seg[peak] = 0.0;
				//seg_heavy[peak] = Atom(0, 0, 0, 0, 0);
			}else{
				done[j] = 2;
				seg_light[peak] = Atom(nl - pos, n, nH, nN, nC);
				for (int u = l; u < r; u++){
					//cout << u << ' ';
					ant.mod_seg[u].done[j] = 1;
					ant.mod_seg[u].seg_num++;
					ant.mod_seg[u].seg_light[peak] = Atom(nl - pos, n, nH, nN, nC);
					ant.mod_seg[u].seg_heavy[peak] = seg_heavy[peak];
					ant.mod_seg[u].mass_seg[peak] = mass_seg[peak];
				}
			}
		}

	}
}

void ModifiedChains::chain_process(){
	string chain = ant.lchain;

	int n = chain.size();
	long double mass = MH2O + MCYST;
	int cyst = -1, pos = 0;

	while (cyst == -1 && pos < n){
		char ch = chain[pos];
		mass += am_wght[ch];
		if (ch == 'C') cyst++;
		pos++;
	}
	if (cyst == -1) return;
	pos--; mass -=am_wght['C']; cyst--;

	for(int i = pos; i < n; i++){

		if (!(i % 50)) cout << i << endl;

		char ch = chain[i];
		if (ch == 'C'){
			cyst++;
		}
		mass += am_wght[ch];
		
		forn(ic, cyst / 2 + 1){

			mass += MCYST * ic;

			heavy_search(mass, 0, 0, ic, i);
			heavy_search(mass - MH2O, 1, 0, ic, i);
			heavy_search(mass - MNH3, 0, 1, ic, i);
			heavy_search(mass - 2 * MH2O, 2, 0, ic, i);
			heavy_search(mass - MH2O - MNH3, 1, 1, ic, i);

			mass -= MCYST * ic;
		}
	}

	forn(i, done.size())
		if (done[i]) seg_num++;

}