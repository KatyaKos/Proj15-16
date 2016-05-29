#include "../headers/headers.h"

using namespace std;

ifstream fin;
ofstream fout;

unordered_map<char, long double> am_wght;
vector<long double> peaks;

bool comp(const ModifiedChains& a, const ModifiedChains& b){
	if (a.nhcyst == 4) return 1;
	if (b.nhcyst == 4) return 0;
	return (a.seg_num > b.seg_num);
}

int main(){

	Reader rd;
	rd.Read();

	/*long double mass = -2*MH2O+MCYST;
	string s = "CSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNV";
	for (char c : s) mass += am_wght[c];
	cout << mass << endl;*/

	string lchain, hchain;
	fin.open(CHAINS_FILE);
	getline(fin, lchain);
	getline(fin, lchain);
	getline(fin, hchain);
	getline(fin, hchain);
	fin.close();

	lchain.pop_back();
	hchain.pop_back();

	Antibody ant(lchain, hchain);
	ant.Calculate();

	Printer pr(ant);
	//pr.Annotate(comp);
	pr.Pict_Annotate(comp);
	//pr.Segments_Cover(comp);
	//pr.Modified_Annotate(comp);

	return 0;
}