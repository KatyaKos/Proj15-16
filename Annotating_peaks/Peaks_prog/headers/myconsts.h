#ifndef _MYCONSTS_H_
#define _MYCONSTS_H_

#define forn(i, n) for(int i = 0; i < n; i++)

const long double MNH3 = 17.03052;
const long double MH2O = 18.01528;
const long double MCYST = -2.01565;

extern long double PREC;

const unordered_map<char, long double> am_wght = {
	{'G', 57.02146},
	{'A', 71.03711},
	{'S', 87.03203},
	{'P', 97.05276},
	{'V', 99.06841},
	{'T', 101.04768},
	{'C', 103.00919},
	{'I', 113.08406},
	{'L', 113.08406},
	{'N', 114.04293},
	{'D', 115.02694},
	{'Q', 128.05858},
	{'K', 128.09496},
	{'E', 129.04259},
	{'M', 131.04049},
	{'H', 137.05891},
	{'F', 147.06841},
	{'R', 156.10111},
	{'Y', 163.06333},
	{'W', 186.07931}	
};

#endif
