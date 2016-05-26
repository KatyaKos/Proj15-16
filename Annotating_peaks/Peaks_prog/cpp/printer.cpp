#include "../headers/headers.h"

using namespace std;

inline string Printer::lost_atoms(int h, int n, int c){
	string lost = "";
	if (h == 1) lost += "-H2O";
	else if (h > 1) lost += "-" + to_string(h) + "H2O";
	if (n == 1) lost += "-NH3";
	else if (n > 1) lost += "-" + to_string(n) + "NH3";
	if (c == 1) lost += "-(S-S)";
	else if (c > 1) lost += "-" + to_string(c) + "(S-S)";

	return lost;
}

inline void Printer::print_seg_peak(const Atom& at, int n){
	int i = at.seg.first, j = at.seg.second;
	if (i){
		string lost = lost_atoms(at.numH2O, at.numNH3, at.numCyst);			

		if (i == 1) fout << setw(36) << right << "b" + to_string(j) + lost;
		else if (j == n) fout << setw(36) << right << "y" + to_string(n - i + 1) + lost;
		else fout << setw(36) << right << "i(" + to_string(i) + '-' + to_string(j) + ')' + lost;

	}else fout << setw(36) << right << "None";
}

inline void Printer::print_mod_peak(const Atom& atl, const Atom& ath){
	int nl = ant.lchain.size();
	string lost = lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ath.numCyst + atl.numCyst);

	fout << setw(36) << right << "ly" + to_string(nl - atl.seg.first + 1) + "+S-S+hi(" + to_string(ath.seg.first) + '-' + to_string(ath.seg.second) + ')' + lost;
			
}

inline void Printer::print_pict_peak(const Atom& at, const string& chain){
	int i = at.seg.first, j = at.seg.second;
	if (i){
		fout << chain.substr(i - 1, j - i + 1);
	}
	
}

void Printer::connect(bool (*comp)(const ModifiedChains&, const ModifiedChains&)){
	sort(ant.mod_seg.begin(), ant.mod_seg.end(), comp);
	int nm = ant.mod_seg.size();


	forn(i, peaks.size()){
		long double peak = peaks[i];

		forn(j, nm){
			if(ant.mod_seg[j].done[i]){
				mod_pos[i] = j;
				ant.mod_num++;
				break;
			}
		}
	}
}

void Printer::where_is_cyst(int ci){
	string lchain = ant.lchain, hchain = ant.hchain;
	int i = ant.posC_heavy[ci], psz = ant.posC_heavy.size(), lsz = lchain.size(), hsz = hchain.size();

	int n = 0;
	forn(u, mod_pos.size()){
		int j = mod_pos[u];
		if (j == -1) continue;

		long double peak = peaks[u];
		Atom ath = ant.mod_seg[j].seg_heavy[peak], atl = ant.mod_seg[j].seg_light[peak];
		int lh = ath.seg.first, rh = ath.seg.second, ll = atl.seg.first;

		if (lh <= i && rh >= i){
			n++;
			string lost = lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ath.numCyst + atl.numCyst);

			fout << peak << ":  ly" + to_string(ll) + "+S-S+hi(" + to_string(lh) + '-' + to_string(rh) + ')' + lost << endl;
			fout << "LC:  ";
			int lcyst = -1, hcyst = -1, ncyst = ath.numCyst + atl.numCyst;
			for (int i = lsz - ll; i < lsz; i++){
				fout << lchain[i];
				if (lchain[i] == 'C') lcyst++;
			}
			fout << endl << "HC:  ";
			for (int i = lh - 1; i <= rh - 1; i++){
				fout << hchain[i];
				if (hchain[i] == 'C') hcyst++;
			}

			fout << endl << "LC mods:  ";
			if (!hcyst && ncyst){
				if (ncyst == 1) fout << "-(S-S)";
				else fout << '-' << ncyst << "(S-S)";
				ncyst = 0; lcyst = 0;
			}else if (lcyst && ncyst && ncyst * 2 == hcyst + lcyst){
				if (lcyst == 2) fout << "-(S-S)";
				else fout << '-' << lcyst / 2 << "(S-S)";
				ncyst = hcyst / 2; lcyst = 0;
			}
			fout << endl << "HC mods:  ";
			if (!lcyst && ncyst){
				if (ncyst == 1) fout << "-(S-S)";
				else fout << '-' << ncyst << "(S-S)";
				ncyst = 0;
			}
			fout << endl << "Other mods:  " << lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ncyst);
			fout << endl << endl;
		}
	}
	fout << endl << "IN TOTAL: " << n << endl;
}

void Printer::lonely_cyst(int ci){

	int i = ant.posC_heavy[ci], psz = ant.posC_heavy.size();
	int n = 0;

	forn(u, mod_pos.size()){
		int k = mod_pos[u];
		if (k == -1) continue;

		long double peak = peaks[u];
		Atom ath = ant.mod_seg[k].seg_heavy[peak], atl = ant.mod_seg[k].seg_light[peak];
		int flag = 0, lpos = ath.seg.first, rpos = ath.seg.second;

		if (lpos <= i && rpos >= i){
			for (int j = 0; j < psz; j++)
				if (ant.posC_heavy[j] <= rpos && ant.posC_heavy[j] >= lpos && j != ci){
					flag = 1;
					break;
				}

			if (flag) continue;
			string lost = lost_atoms(ath.numH2O + atl.numH2O, ath.numNH3 + atl.numNH3, ath.numCyst + atl.numCyst);

			fout << setw(36) << right << "ly" + to_string(atl.seg.first) + "+S-S+hi(" + to_string(ath.seg.first) + '-' + to_string(ath.seg.second) + ')' + lost << endl;
			n++;
		}
	}
	fout << endl << "IN TOTAL: " << n << endl;

}

pair<int, int> Printer::LRtest(int ci){
	int i = ant.posC_heavy[ci], psz = ant.posC_heavy.size();
	int nseg = 0, n = 0, no = 0;

	forn(u, mod_pos.size()){
		int k = mod_pos[u];
		if (k == -1) continue;

		long double peak = peaks[u];
		Atom ath = ant.mod_seg[k].seg_heavy[peak], atl = ant.mod_seg[k].seg_light[peak];
		int flag1 = 1, flag2 = 1;
		int lpos = ath.seg.first, rpos = ath.seg.second;

		if (lpos <= i && rpos >= i){
			n++;
			for (int j = 0; j < psz; j++){
				int jpos = ant.posC_heavy[j];
				if (jpos < i && jpos >= lpos) flag1 = 0;	
				if (j != ci && jpos >= lpos && jpos <= rpos) flag2 = 0;

			}
			if (flag1) nseg++;
			if (flag2) no++;
		}
	}

	if (ci == 0 || ci == psz - 1) return make_pair(n, no);
	return make_pair(nseg, -1);
}


void Printer::print_cyst_group(int nC){
	fout << endl << "PEAK " << nC + 1 << endl;
	fout << "------------------------------------------" << endl;

	ModifiedChains mc = ant.mod_seg[nC];
	um_lda seg_light = mc.seg_light;
	um_lda seg_heavy = mc.seg_heavy;
	int nsz = peaks.size();
	vector<long double> None, OneH, OneN, HandN, TwoH, Other;

	forn(i, nsz){
		if (mc.done[i]){
			long double peak = peaks[i];
			Atom ath = seg_heavy[peak], atl = seg_light[peak];
			int nH = ath.numH2O + atl.numH2O, nN = ath.numNH3 + atl.numNH3;
			if (!nH && !nN) None.push_back(peak);
			else if (nH == 1 && !nN) OneH.push_back(peak);
			else if (!nH && nN == 1) OneN.push_back(peak);
			else if (nH == 1 && nN == 1) HandN.push_back(peak);
			else if (nH == 2 && !nN) TwoH.push_back(peak);
			else Other.push_back(peak);
		}
	}

	long double peak;
	forn(i, None.size()){
		peak = None[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	forn(i, OneH.size()){
		peak = OneH[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	forn(i, OneN.size()){
		peak = OneN[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	forn(i, HandN.size()){
		peak = HandN[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	forn(i, TwoH.size()){
		peak = TwoH[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	forn(i, Other.size()){
		peak = Other[i];
		fout << peak << ":  "; print_mod_peak(seg_light[peak], seg_heavy[peak]); fout << endl;
	}
	
	fout << "------------------------------------------" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Printer::Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&)){

	fout.open(ANNOT_FILE);

	fout << "\nFirst column -- peaks, second -- light chain segment, third -- heavy chain segment, fourth -- modified segment\n\n";

	int nl = ant.lchain.size(), nh = ant.hchain.size(), nm = ant.mod_seg.size();
	connect(comp);

	forn(i, peaks.size()){
		long double peak = peaks[i];
		fout.precision(10);
		fout << setw(11) << right << peak;

		Atom ath_s = ant.hseg.seg[peak], atl_s = ant.lseg.seg[peak];

		print_seg_peak(atl_s, nl);
		print_seg_peak(ath_s, nh);

		int j = mod_pos[i];
		if(j != -1){
			Atom ath_m = ant.mod_seg[j].seg_heavy[peak], atl_m = ant.mod_seg[j].seg_light[peak];
			print_mod_peak(atl_m, ath_m);
		} else fout << setw(36) << right << "None";
		fout << endl;

		fout << setw(11) << right << "DELTA MASS=";
		fout << setw(36) << right << peak - ant.lseg.mass_seg[peak];
		fout << setw(36) << right << peak - ant.hseg.mass_seg[peak];
		if (j != -1) fout << setw(36) << right << peak - ant.mod_seg[j].mass_seg[peak];
		else fout << setw(36) << right << peak;
		fout << endl << endl;

	}

	fout << "LIGHT_SEGMENTS_FOUND=" << ant.light_num << endl;
	fout << "HEAVY_SEGMENTS_FOUND=" << ant.heavy_num << endl;
	fout << "MODIFIED_SEGMENTS_FOUND=" << ant.mod_num << endl;

	fout.close();

}

void Printer::Pict_Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&)){

	fout.open(PICT_FILE);

	fout << "\nOne block -- one peak. First raw in block -- peak, second-third -- light chain segment, fourth-fifth -- heavy chain segment, sixth-seventh -- modified segment\n\n";

	int nl = ant.lchain.size(), nh = ant.hchain.size();
	connect(comp);

	forn(i, peaks.size()){
		fout << endl;
		long double peak = peaks[i];
		fout.precision(10);
		fout << "--------------------" << endl;
		fout << "PEAK:  " << peak << endl << endl;

		Atom ath_s = ant.hseg.seg[peak], atl_s = ant.lseg.seg[peak];

		fout << "LIGHT FOUND:     "; print_seg_peak(atl_s, nl); fout << endl;
		print_pict_peak(atl_s, ant.lchain); fout << endl;
		fout << "DELTA MASS = " << peak - ant.lseg.mass_seg[peak] << endl << endl;

		fout << "HEAVY FOUND:     "; print_seg_peak(ath_s, nh); fout << endl;
		print_pict_peak(ath_s, ant.hchain); fout << endl;
		fout << "DELTA MASS = " << peak - ant.hseg.mass_seg[peak] << endl << endl;

		fout << "MODIFIED FOUND:  "; 
		int j = mod_pos[i];
		if(j != -1){

			Atom ath_m = ant.mod_seg[j].seg_heavy[peak], atl_m = ant.mod_seg[j].seg_light[peak];
			print_mod_peak(atl_m, ath_m);
			fout << endl;
			print_pict_peak(atl_m, ant.lchain); fout << "<-->"; print_pict_peak(ath_m, ant.hchain);
		
		}else fout << setw(36) << right << "None" << endl;
		
		fout << endl;
		if (j != -1) fout << "DELTA MASS = " << peak - ant.mod_seg[j].mass_seg[peak] << endl;
		else fout << "DELTA MASS = " << peak << endl;

		fout << "--------------------" << endl << endl;
	}

	fout << "LIGHT_SEGMENTS_FOUND=" << ant.light_num << endl;
	fout << "HEAVY_SEGMENTS_FOUND=" << ant.heavy_num << endl;
	fout << "MODIFIED_SEGMENTS_FOUND=" << ant.mod_num << endl;

	fout.close();

}

void Printer::Segments_Cover(bool (*comp)(const ModifiedChains&, const ModifiedChains&)){
	fout.open(CYS_PROCESS_FILE);

	connect(comp);

	fout << "MODIFIED SEGMENTS THAT CONTAIN FIFTH CYSTEIN:" << endl << endl;
	where_is_cyst(4);
	fout << endl << endl;
	fout << "MODIFIED SEGMENTS THAT CONTAIN ONLY FIFTH CYSTEIN:" << endl << endl;
	lonely_cyst(4);

	int psz = ant.posC_heavy.size();
	fout << endl << endl << endl << "NUMBER OF MODIFIED SEGMENTS THAT COVER CYSTEINS WITH THEIR LEFTEST/RIGHTEST END:" << endl;
	for (int i = 1; i < psz - 1; i++){
		fout << "cystein " << i + 1 << ":  " <<  LRtest(i).first << endl;
	}

	pair<int, int> pf = LRtest(0), pl = LRtest(psz - 1);
	fout << endl << endl << "NUMBER OF MODIFIED SEGMENTS THAT CONTAIN Nth CYSTEIN:" << endl;
	fout << "cystein 1:  " << pf.first << endl << "cystein " << psz << ":  " << pl.first << endl << endl;
	fout << endl << endl << "NUMBER OF MODIFIED SEGMENTS THAT CONTAIN ONLY Nth CYSTEIN:" << endl;
	fout << "cystein 1:  " << pf.second << endl << "cystein " << psz << ":  " << pl.second << endl << endl;

	fout.close();
}

void Printer::Modified_Annotate(bool (*comp)(const ModifiedChains&, const ModifiedChains&)){
	fout.open(MOD_FILE);

	connect(comp);
	int nC = ant.posC_heavy.size();

	print_cyst_group(4);
	forn(i, nC){
		if (i == 4) continue;
		print_cyst_group(i);
	}
	fout.close();
}