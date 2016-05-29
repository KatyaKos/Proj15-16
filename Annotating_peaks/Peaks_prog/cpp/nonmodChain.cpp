#include "../headers/headers.h"

using namespace std;

void NonModChains::chain_process(){
	done.assign(peaks.size(), 0);

	forn(i, chain.size()){
		peak_search(seg, done, chain, 0.0, i);
	}

}