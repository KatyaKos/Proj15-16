#ifndef _METHODS_H_
#define _METHODS_H_

//work with info
void Read_weights();
void Read_peaks();
void annotating(const string& lchain, const string& hchain);
void pict_annotating(const string& lchain, const string& hchain);
void modified_annotating(const string& lchain, const string& hchain);

//parse keys
void chains_modified(string lchain, string hchain);
int chains_normal(const string& chain, um_lda& annot, um_ldld& anmass);
void SegCover(const string& lchain, const string& hchain);

//work with peaks
int check(long double m, int j);

int C_preprocess(const string& chain, vector<int>& place);
int modified_peak_search(vector<int>& done, long double delta, int cyst, const string& lchain, const string& hchain, int pos, int posC);

int peak_search(um_lda& annot, um_ldld& anmass, vector<int>& done, long double delta, const string& chain, int pos);

//tell me more about cysteins in modified peaks
void Where_is_cyst(int ci, const string& lchain, const string& hchain);
void LonelyCyst(int ci);
pair<int, int> LRtest(int ci);

#endif