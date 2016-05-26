#ifndef _METHODS_H_
#define _METHODS_H_

//work with info
void pict_annotating(const string& lchain, const string& hchain);
void modified_annotating(const string& lchain, const string& hchain);

//parse keys
void SegCover(const string& lchain, const string& hchain);

//tell me more about cysteins in modified peaks
void Where_is_cyst(int ci, const string& lchain, const string& hchain);
void LonelyCyst(int ci);
pair<int, int> LRtest(int ci);

#endif