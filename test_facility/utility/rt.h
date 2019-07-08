#include<iostream>
#include<map>
#include<set>
using namespace std;

#define CLS 64
#define DS 8

/* first access time */
std::map<int, int> fat;
/* last access time */
std::map<int, int> lat;
/* reuse time histogram */
std::map<uint64_t, uint64_t> rtTmp;
unsigned long long refT = 0;

std::map<uint64_t, double> RT;

std::map<uint64_t, double> MR;

std::map<uint64_t, std::set<uint64_t> > rtRefSet;
std::map<uint64_t, std::set<uint64_t> > rtArrSet;

void rtTmpAccess(uint64_t addr, uint64_t ref_id, uint64_t array_id) {
	addr = addr * DS / CLS;

	refT++;
	if (lat.find(addr) == lat.end()) {
        //fat[addr] = refT;
        lat[addr] = refT;
    } else {
		/*
        if (rtTmp.find(refT - lat[addr]) == rtTmp.end()) {
            rtTmp[refT - lat[addr]] = 1;
		} else {
            rtTmp[refT - lat[addr]] ++;
		}
		*/
		rtRefSet[ref_id].insert(refT - lat[addr]);
		rtArrSet[array_id].insert(refT - lat[addr]);
        lat[addr] = refT;
    }
	return;	
}

void rtTmpAccess(int addr) {

    addr = addr * DS / CLS;

    refT++;
    if (fat.find(addr) == fat.end()) {
        fat[addr] = refT;
        lat[addr] = refT;
    } else {

        if (rtTmp.find(refT - lat[addr]) == rtTmp.end()) {
            rtTmp[refT - lat[addr]] = 1;
        } else {
            rtTmp[refT - lat[addr]] ++;
        }
        lat[addr] = refT;
    }
    return;
}

void dumpRtTmp() {
    uint64_t cnt = 0;
    for (std::map<uint64_t, uint64_t>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cnt += it->second;
    }

    for (std::map<uint64_t, uint64_t>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cout << it->first << " " << it->second << " " << double(it->second)/cnt  << endl;
    }

	cout << rtTmp.size() << endl;
    //rtTmp.clear();
    //fat.clear();
    //lat.clear();

    return;
}

void RTtoMR_AET() {

	for (std::map<uint64_t, uint64_t>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
		RT[it->first] = (double) it->second;
	}

	std::map<uint64_t, double> P;
	double total_num_RT = 0;
	uint64_t max_RT = 0;

	for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
		total_num_RT += it->second;
		if (max_RT < it->first) {
			max_RT = it->first;
		}
	}

	double accumulate_num_RT = 0;
	for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
		P[it->first] = accumulate_num_RT / total_num_RT;
		//cout << accumulate_num_RT << " " << total_num_RT << endl;
		accumulate_num_RT += it->second;
		//cout << "P " << it->first << " " << P[it->first] << endl;
	}

	P[0] = 1;

	double sum_P = 0;
	uint64_t t = 0;
	uint64_t prev_t = 0;
	
	double MR_pred = -1;


//	std::cout << "here " << std::endl;

	// 20MB
	for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {
//	for (uint64_t c = 0; c <= max_RT; c++) {
		while (sum_P < c && t <= max_RT) {
			if (P.find(t) != P.end()) {
				sum_P += P[t];
				prev_t = t;
			} else {
				sum_P += P[prev_t];
			}
			t++;
		}

		
		if (MR_pred != -1) {
			MR[c] = P[prev_t];
			MR_pred = P[prev_t];
		} else {
			if (MR_pred - P[prev_t] < 0.0001) {
				MR[c] = P[prev_t];
            	MR_pred = P[prev_t];
			}
		}		


		//MR[c] = P[prev_t];

	}
	
	return;
}


void dumpMR() {
	

	cout << "miss ratio" << endl;	
	std::map<uint64_t, double>::iterator it1 = MR.begin();
	std::map<uint64_t, double>::iterator it2 = MR.begin();

	while(it1 != MR.end()) {
		while(1) {
			std::map<uint64_t, double>::iterator it3 = it2;
			++it3;
			if (it3 == MR.end()) {
				break;
			}
			//if (it3->second == it2->second) {
			if (it1->second - it3->second < 0.00001) {
				++it2;
			} else {
				break;
			}
		}
		cout << it1->first << " " << it1->second << endl;
		if (it1 != it2) {
			cout << it2->first << " " << it2->second << endl;
		}
		it1 = ++it2;
		it2 = it1;
	}

	return;
}

void dumpSetSize() {
	cout << "#different RT value Per reference:" << endl;
	uint64_t cnt = 0;
	double size = 0;
	for (std::map<uint64_t, std::set<uint64_t> >::iterator it = rtRefSet.begin(), eit = rtRefSet.end(); it != eit; ++it) {
		cout << it->second.size() << endl;
		size += it->second.size();
		cnt ++;
	}
	cout << "Average " << size / cnt << endl << endl;
	
	cout << "#different RT value Per array:" << endl;
	cnt = 0;
	size = 0;
	for (std::map<uint64_t, std::set<uint64_t> >::iterator it = rtArrSet.begin(), eit = rtArrSet.end(); it != eit; ++it) {
		cout << it->second.size() << endl;
        size += it->second.size();
        cnt ++;
	}
	cout << "Average " << size / cnt << endl << endl;
}

