#ifndef _ANTIBODY_H_
#define _ANTIBODY_H_

using namespace std;

struct Atom{
	int numH2O, numNH3, numCyst;
	pair <int, int> seg;
	Atom() = default;
	Atom(int i, int j, int h, int n, int c): seg(make_pair(i, j)), numCyst(c), numNH3(n), numH2O(h) {}
};

typedef unordered_map<long double, Atom> um_lda;
typedef unordered_map<long double, long double> um_ldld;

class Antibody;

class Chains{
protected:
	const Antibody& ant;
	int check(long double m, int j);
public:
	int seg_num;
	vector<int> done;

	um_ldld mass_seg;

	Chains(const Antibody& ant): ant(ant), seg_num(0) {}
};

class NonModChains: public Chains{
private:
	string chain;
	void peak_search(long double delta, int pos);
public:
	um_lda seg;

	NonModChains() = default;
	NonModChains(const Antibody& ant, string chain): Chains(ant), chain(chain) {}
	NonModChains(const NonModChains& c): Chains(c.ant){
		done = c.done;
		chain = c.chain;
		seg = c.seg;
		seg_num = c.seg_num;
		mass_seg = c.mass_seg;
	}

	void chain_process();
};

class ModifiedChains: public Chains{
private:
	void peak_search(long double delta, int cyst, int pos, int nC);
	int light_peak_search(long double delta, int pos, int j);
public:
	int nhcyst;

	um_lda seg_light;
	um_lda seg_heavy;

	ModifiedChains() = default;
	ModifiedChains(const Antibody& ant, int i): Chains(ant), nhcyst(i) {}
	ModifiedChains& operator=(const ModifiedChains& c){
		Chains(c.ant);
		done = c.done;
		seg_num = c.seg_num, nhcyst = c.nhcyst;
		seg_heavy = c.seg_heavy, seg_light = c.seg_light;
		mass_seg = c.mass_seg;
	}

	void chain_process(int ii);
};

class Antibody{
private:
	void C_preprocess();
public:
	string lchain, hchain;
	int light_num, heavy_num, mod_num;
	vector<int> posC_heavy, posC_light;

	vector<ModifiedChains> mod_seg;
	NonModChains lseg;
	NonModChains hseg;

	Antibody() = default;
	Antibody(string lch, string hch): light_num(0), heavy_num(0), mod_num(0), lchain(lch), hchain(hch), lseg(*this, lch), hseg(*this, hch) {}
	Antibody(const Antibody& ant): lchain(ant.lchain), hchain(ant.hchain), lseg(ant.lseg), hseg(ant.hseg) {
		posC_heavy = ant.posC_heavy, posC_light = ant.posC_light;
		mod_seg = ant.mod_seg;
		light_num = ant.light_num, heavy_num = ant.heavy_num, mod_num = ant.mod_num;
	}

	void Calculate();
};

#endif