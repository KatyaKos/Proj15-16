#include "../headers/headers.h"

using namespace std;

int Chains::check(long double m, int j){
	long double peak = peaks[j];

	if (m <= peak * 1.00001 && m >= peak * 0.99999) return 1;
	return 0;
}