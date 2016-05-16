#ifndef _MYCONSTS_H_
#define _MYCONSTS_H_

#define AMINO_FILE "/home/student/BioinfProj/Proj15-16/Annotating_peaks/aminos.in"
#define SPECTRUM_FILE "/home/student/BioinfProj/Proj15-16/Annotating_peaks/mass_spectrum.msalign"
#define CHAINS_FILE "/home/student/BioinfProj/Proj15-16/Annotating_peaks/trastu-fab.fasta"
#define RESULT_FILE "/home/student/BioinfProj/Proj15-16/Annotating_peaks/Peaks_prog/annotated_peaks.txt"
#define CYS_PROCESS_FILE "/home/student/BioinfProj/Proj15-16/Annotating_peaks/Peaks_prog/cys_process.txt"

const int BUFFN = 50;
const long double MNH3 = 17.03052;
const long double MH2O = 18.01528;
const long double MCYST = -2.01565;
const long double PREC = 0.00001;

//long double PRECURSOR_MZ;                   //это нужно?
//int PRECURSOR_CHARGE;
//long double PRECURSOR_MASS;

#endif