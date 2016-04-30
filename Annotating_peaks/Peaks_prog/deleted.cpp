inline void check_maps(){

	fout.open("/home/student/Practice/Task_2/all_peaks.txt");

	fout << "LIGHT_PREFIX_FOUND=" << light_pref_num << endl;
	fout << "LIGHT_SUFFIX_FOUND=" << light_suf_num << endl;
	fout << "HEAVY_PREFIX_FOUND=" << heavy_pref_num << endl;
	fout << "HEAVY_SUFFIX_FOUND=" << heavy_suf_num << endl;
	fout << "HEAVY_SEGMENTS_FOUND=" << light_seg_num << endl;
	fout << "HEAVY_SEGMENTS_FOUND=" << heavy_seg_num << endl;
	fout << "\nFirst column -- peaks, second -- light prefs/sufs/segmets, third -- heavy prefs/sufs/segmets.\n\n";

	forn(i, peaks.size()){
		long double peak = peaks[i];
		fout.precision(10);
		fout.width(11);
		fout << peak;

		if (pref_light[peak].first.second)
			fout << setw(40) << right << "b" + to_string(pref_light[peak].first.second) + pref_light[peak].second;
		
		else if (suf_light[peak].first.second)
			fout << setw(40) << right << "y" + to_string(suf_light[peak].first.second) + suf_light[peak].second;
		
		else if (seg_light[peak].first.second)
			fout << setw(40) << right << "i(" + to_string(seg_light[peak].first.first) + '-' + to_string(seg_light[peak].first.second) + ')' + seg_light[peak].second;

		else fout << setw(40) << right << "None";
		fout << "        ";

		if (pref_heavy[peak].first.second)
			fout << setw(40) << right << "b" + to_string(pref_heavy[peak].first.second) + pref_heavy[peak].second;

		else if (suf_heavy[peak].first.second)
			fout << setw(40) << right << "y" + to_string(suf_heavy[peak].first.second) + suf_heavy[peak].second;

		else if(seg_heavy[peak].first.second)
			fout << setw(40) << right << "i(" + to_string(seg_heavy[peak].first.first) + '-' + to_string(seg_heavy[peak].first.second) + ')' + seg_heavy[peak].second;

		else fout << setw(40) << right << "None";
		fout << endl;
	}
	fout.close();

}

int main(){

	Read_weights();
	Read_peaks();

	chains_process();

	string s = "TNGYTR";
	long double tmp = - MNH3;
	forn(i, s.size()){
		tmp += am_wght[s[i]];
	}
	cout << tmp << '\n';

	check_maps();

	return 0;
}