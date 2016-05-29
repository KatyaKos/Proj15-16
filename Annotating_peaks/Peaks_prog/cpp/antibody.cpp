#include "../headers/headers.h"

using namespace std;

void Antibody::C_preprocess(){
	forn(i, lchain.size())
		if (lchain[i] == 'C')
			posC_light.push_back(i);

	forn(i, hchain.size())
		if (hchain[i] == 'C')
			posC_heavy.push_back(i);
		
}

void Antibody::Calculate(){
	lseg.chain_process();
	hseg.chain_process();

	light_num = lseg.seg_num, heavy_num = hseg.seg_num;

	C_preprocess();
	int nhC = posC_heavy.size();
	if (!nhC) return;
	forn(ii, nhC)
		mod_seg.push_back(ModChains(*this, ii));

	reverse(lchain.begin(), lchain.end());
	unite.chain_process();
	reverse(lchain.begin(), lchain.end());

}