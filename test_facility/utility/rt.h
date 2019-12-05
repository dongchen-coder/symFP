#include<iostream>
#include<map>
#include<set>
#include <cmath>
using namespace std;

#define CLS 64
#define DS 8
#ifndef BIN_SIZE
#    define BIN_SIZE   64
#endif

/* first access time */
std::map<int, int> fat;
/* last access time */
std::map<int, int> lat;
/* reuse time histogram */
std::map<uint64_t, uint64_t> rtTmp;
unsigned long long refT = 0;

std::map<uint64_t, double> RT;

std::map<uint64_t, double> MR;

map<uint64_t, map<uint64_t, uint64_t>* > RI;

map<uint64_t, uint64_t> refAccessCnt;
map<uint64_t, uint64_t> srcRef;

//std::map<uint64_t, std::set<uint64_t> > rtRefSet;
//std::map<uint64_t, std::set<uint64_t> > rtArrSet;
void rtHistoCal( uint64_t rt, uint64_t val ) {
    if ( val <= 0) {
;        return;
    }
    if (rtTmp.find(rt) == rtTmp.end()) { 
        rtTmp[rt] = val;
    } else {
        rtTmp[rt] += val;
    }
    return;
}
void subBlkRT(uint64_t rt) {
    int msb = 0;
    int tmp_rt = rt;
    while(tmp_rt != 0) {
        tmp_rt = tmp_rt / 2;
        ++msb;
    }
    if (msb >= BIN_SIZE) {
        uint64_t diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;
        for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {
            if (rt < b) {
                rtHistoCal((uint64_t)(b - diff), 1);
                break;
            }
        }
    }
    else {
        rtHistoCal((uint64_t)pow(2, msb-1), 1);
    }
    return;
}

void rtTmpAccess(uint64_t addr, uint64_t ref_id, uint64_t array_id) {
	addr = addr * DS / CLS;
	refT++;

	if (lat.find(addr) != lat.end()) {
        uint64_t ri = refT - lat[addr];
        uint64_t sourceRef = srcRef[addr];
        if (RI.find(sourceRef) != RI.end()) {
            if ((*RI[sourceRef]).find(ri) != (*RI[sourceRef]).end()) {
                (*RI[sourceRef])[ri] ++;
            } else {
                (*RI[sourceRef])[ri] = 1;
            }
        } else {
            RI[sourceRef] = new map<uint64_t, uint64_t>;
            (*RI[sourceRef])[ri] = 1;
        }
    }
    lat[addr] = refT;
    srcRef[addr] = ref_id;

	if (refAccessCnt.find(ref_id) != refAccessCnt.end()) {
		refAccessCnt[ref_id] ++;
	} else {
		refAccessCnt[ref_id] = 1;
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
    	subBlkRT(refT - lat[addr]);
        lat[addr] = refT;
    }
    return;
}

void dumpRtTmp() {
	cout << "Number of Cache Line: " << fat.size() * DS / CLS << endl;
    uint64_t cnt = 0;
    for (std::map<uint64_t, uint64_t>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cnt += it->second;
    }
    cout << "Start to dump reuse time histogram" << endl;
    for (std::map<uint64_t, uint64_t>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cout << it->first << ", " << it->second << ", " << double(it->second)/cnt  << endl;
    }

	// cout << rtTmp.size() << endl;

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

void dumpRI() {
	uint64_t total_number_of_ri = 0;
	uint64_t total_number_of_ri_cnt = 0;
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		std::set<uint64_t> riset;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			cout << "Ref " << ref_it->first << " RI " << ri_it->first << " CNT " << ri_it->second << endl;
			riset.insert(ri_it->first);
			total_number_of_ri_cnt += ri_it->second;
		}
		cout << "Ref " << ref_it->first << " RISETSIZE " << riset.size() << endl;
		total_number_of_ri += riset.size();
	}
	cout << "Average RISETSIZE for each reference " << double(total_number_of_ri) / RI.size() << endl;
	
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            cout << "Ref " << ref_it->first << " RI " << ri_it->first << " DIST " << double(ri_it->second) / total_number_of_ri_cnt << endl;
        }
    }

	cout << "Total number of accesses " << refT << endl;
	for (map<uint64_t, uint64_t>::iterator it = refAccessCnt.begin(), eit = refAccessCnt.end(); it != eit; ++it) {
		cout << "Ref " << it->first << " ACCESSCNT " << it->second << endl;
	}
}

