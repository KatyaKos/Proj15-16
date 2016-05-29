#include "../headers/headers.h"

using namespace std;

ifstream fin;
ofstream fout;

long double PREC;
vector<long double> peaks;

void help(){
	cerr << "ERROR!" << endl;
    cerr << "Usage: ./prog -f *.fasta -m *.msalign -p PRECISION" << endl;
}

void check_args_num(int argc){
	if (argc != 1+6) throw "There must be three keys.";
}

void check_args(int flag){
	if (flag) throw "There must be keys -f, -o, -c/-u.";
}

void check_file(string name){
	if (ifstream(name, ios::in) == NULL) throw "File problem";
}



bool comp(const ModifiedChains& a, const ModifiedChains& b){
	if (a.nhcyst == 4) return 1;
	if (b.nhcyst == 4) return 0;
	return (a.seg_num > b.seg_num);
}


int main(int argc, char* argv[]){

	try{
		check_args_num(argc);
	} catch (const char* s){
		help();
		return 1;
	}

	string spect, anti;
	int flag = 1;

	for (int i = 1; i < 6; i++)
		if (!strcmp(argv[i], "-m"))
			spect = argv[i + 1], flag = 0;
		
	try{
		check_args(flag);
	} catch (const char* s){
		help();
		return 1;
	}
	flag = 1;

	for (int i = 1; i < 6; i++)
		if (!strcmp(argv[i], "-f"))
			anti = argv[i + 1], flag = 0;
		
	try{
		check_args(flag);
	} catch (const char* s){
		help();
		return 1;
	}
	flag = 1;

	for (int i = 1; i < 6; i++)
		if (!strcmp(argv[i], "-p")){
			char* tmparg = argv[i + 1];
			char* ptEnd;
			PREC = strtold(tmparg, &ptEnd);
			flag = 0;
		}
		
	try{
		check_args(flag);
	} catch (const char* s){
		help();
		return 1;
	}

	try{
		check_file(spect);
		check_file(anti);
	} catch (const char* s) {
		cerr << "ERROR!\nI/O error" << endl;
		return 1;
	}


	Spectrum sp(spect, anti);
	sp.Read();

	/*long double mass = -2*MH2O+MCYST;
	string s = "CSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNV";
	for (char c : s) mass += am_wght.find(c)->second;
	cout << mass << endl;*/

	string lchain, hchain;
	fin.open(anti);
	getline(fin, lchain);
	getline(fin, lchain);
	getline(fin, hchain);
	getline(fin, hchain);
	fin.close();

	lchain.pop_back();
	hchain.pop_back();

	Antibody ant(lchain, hchain);
	ant.Calculate();

	sp.set_antibody(ant);
	sp.Annotate(comp);
	sp.Pict_Annotate(comp);
	sp.Segments_Cover(comp);
	sp.Modified_Annotate(comp);

	return 0;
}