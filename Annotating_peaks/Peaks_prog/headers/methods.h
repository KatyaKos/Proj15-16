#ifndef _METHODS_H_
#define _METHODS_H_

//work with info
void Read_weights();
void Read_peaks();
void check_maps(const string& lchain, const string& hchain);

//parse keys
void chains_modified(string lchain, string hchain);
int chains_normal(string chain, um_lda& annot, um_ldld& anmass);
void SegCover();

//work with peaks
int check(long double m, int j);

int C_preprocess(string chain, vector<int>& place);
int modified_peak_search(vector<int>& done, long double delta, int cyst, string lchain, string hchain, int pos, int posC);

int peak_search(um_lda& annot, um_ldld& anmass, vector<int>& done, long double delta, string chain, int pos);

//tell me more about cysteins in modified peaks
void LonelyCyst(int i);
pair<int, int> LRtest(int ci);

#endif