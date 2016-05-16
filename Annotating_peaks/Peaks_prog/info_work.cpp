#include "headers.h"

using namespace std;

void Read_weights(){
	fin.open(AMINO_FILE);
	char name, tmp;
	long double w;
	while (fin >> name){
		fin >> tmp >> w >> tmp;
		am_wght[name] = w;
	}
	fin.close();
}

void Read_peaks(){
	fin.open(SPECTRUM_FILE);
	string buff;
	forn(i, 6){
		getline(fin, buff);
	}
	long double d;

	while (fin >> d){
		peaks.push_back(d);
		fin >> d >> d;                     //эти два запоминать куда-то?
	}

	sort(peaks.begin(), peaks.end());
	fin.close();
}


string Printer::lost_atoms(int h, int n, int c){
	string lost = "";
	if (h == 1) lost += "-H2O";
	else if (h > 1) lost += " -" + to_string(h) + "H2O";
	if (n == 1) lost += "-NH3";
	else if (n > 1) lost += " -" + to_string(n) + "NH3";
	if (c == 1) lost += "+(S-S)";
	else if (c > 1) lost += "-" + to_string(c) + "(S-S)";

	return lost;
}

void Printer::print_seg_peak(Atom at, int n){
	int i = at.seg.first, j = at.seg.second;
	if (i){
		string lost = lost_atoms(at.numH2O, at.numNH3, at.numCyst);			

		if (i == 1) fout << setw(36) << right << "b" + to_string(j) + lost;
		else if (j == n) fout << setw(36) << right << "y" + to_string(i) + lost;
		else fout << setw(36) << right << "i(" + to_string(i) + '-' + to_string(j) + ')' + lost;

	}else fout << setw(36) << right << "None";
}


void check_maps(const string& lchain, const string& hchain){

	fout.open(RESULT_FILE);

	fout << "LIGHT_SEGMENTS_FOUND=" << seg_num_light << endl;
	fout << "HEAVY_SEGMENTS_FOUND=" << seg_num_heavy << endl;
	fout << "MODIFIED_SEGMENTS_FOUND=" << mod_num << endl;
	fout << "\nFirst column -- peaks, second -- light chain segment, third -- heavy chain segment, fourth -- modified segment\n\n";

	int nl = lchain.size(), nh = hchain.size();
	Printer pr;

	forn(i, peaks.size()){
		long double peak = peaks[i];
		fout.precision(10);
		fout << setw(11) << right << peak;

		Atom ath_m = mod_heavy[peak], atl_m = mod_light[peak];
		Atom ath_s = seg_heavy[peak], atl_s = seg_light[peak];

		pr.print_seg_peak(atl_s, nl);
		pr.print_seg_peak(ath_s, nh);

		if(ath_m.seg.first){
			string lost = pr.lost_atoms(ath_m.numH2O + atl_m.numH2O, ath_m.numNH3 + atl_m.numNH3, ath_m.numCyst + atl_m.numCyst);

			fout << setw(36) << right << "ly" + to_string(atl_m.seg.first) + "+S-S+hi(" + to_string(ath_m.seg.first) + '-' + to_string(ath_m.seg.second) + ')' + lost;
		}else fout << setw(36) << right << "None";
		fout << endl;

		fout << setw(11) << right << "DELTA MASS=";
		fout << setw(36) << right << peak - mass_seg_light[peak];
		fout << setw(36) << right << peak - mass_seg_heavy[peak];
		fout << setw(36) << right << peak - mass_mod[peak];
		fout << endl << endl;

	}

	fout.close();

}