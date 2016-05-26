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

	C_preprocess();

	int nhC = posC_heavy.size();
	forn(ii, nhC)
		mod_seg.push_back(ModifiedChains(*this, ii));

	reverse(lchain.begin(), lchain.end());
	forn(ii, nhC){
		mod_seg[ii].chain_process(ii);
	}
	reverse(lchain.begin(), lchain.end());

	light_num = lseg.seg_num, heavy_num = hseg.seg_num;

}