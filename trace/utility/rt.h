#include<iostream>
#include<map>
using namespace std;

#define CLS 32
#define DS 8

/* first access time */
std::map<int, int> fat;
/* last access time */
std::map<int, int> lat;
/* reuse time histogram */
std::map<int, int> rtTmp;
unsigned long long refT = 0;


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
    for (std::map<int, int>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cnt += it->second;
    }

    for (std::map<int, int>::iterator it = rtTmp.begin(), eit = rtTmp.end(); it != eit; ++it) {
        cout << it->first << " " << it->second << " " << double(it->second)/cnt  << endl;
    }

    rtTmp.clear();
    fat.clear();
    lat.clear();

    return;
}


