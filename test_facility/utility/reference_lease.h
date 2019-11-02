#include <iostream>
#include <set>
#include <map>
#include <limits>
using namespace std;

#define CLS 64
#define DS 8

uint64_t refT = 0;

map<uint64_t, uint64_t> lat;

map<uint64_t, map<uint64_t, uint64_t>* > RI;

map<uint64_t, map<uint64_t, double>* > hits;
map<uint64_t, map<uint64_t, double>* > costs;

//map<uint64_t, double> sampledCnt;
//map<uint64_t, double> accessRatio;

map<uint64_t, uint64_t> Lease;
map<uint64_t, uint64_t> srcRef;

std::map<uint64_t, map<uint64_t, bool>* > refAddrReuseMark;

// get reuse interval distribution
void rtTmpAccess(uint64_t addr, uint64_t ref_id, uint64_t arr_id) {
	refT++;
	addr = addr * DS / CLS;
	
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

	if (refAddrReuseMark.find(ref_id) != refAddrReuseMark.end()) {
		if ((*refAddrReuseMark[ref_id]).find(addr) != (*refAddrReuseMark[ref_id]).end()) {
			(*refAddrReuseMark[ref_id])[addr] = true;
		} else {
			(*refAddrReuseMark[ref_id])[addr] = false;
		}
	} else {
		refAddrReuseMark[ref_id] = new map<uint64_t, bool>;
		(*refAddrReuseMark[ref_id])[addr] = false;
	}

	// Init leases to all references to be 0
	if (Lease.find(ref_id) == Lease.end()) {
		Lease[ref_id] = 0;
	}

	return;
}

void RIwithInfinite() {
	for (map<uint64_t, map<uint64_t, bool>* >::iterator ref_it = refAddrReuseMark.begin(), ref_eit = refAddrReuseMark.end(); ref_it != ref_eit; ++ref_it) {
		uint64_t cnt = 0;
		for(map<uint64_t, bool>::iterator addr_it = (*ref_it->second).begin(), addr_eit = (*ref_it->second).end(); addr_it != addr_eit; ++addr_it) {
			if (addr_it->second == false) {
				cnt++;
			}
		}
		(*RI[ref_it->first])[std::numeric_limits<uint64_t>::max()] = cnt;
	}
	return;
}


/*
// calculate access ratio
void accessRatioCal() {
	double total_access_cnt = 0;
	
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			total_access_cnt += ri_it->second;
		}
	}
	
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		double ref_access_cnt = 0;
		for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			ref_access_cnt += ri_it->second;	
		}
		sampledCnt[ref_it->first] = ref_access_cnt;
		accessRatio[ref_it->first] = ref_access_cnt / total_access_cnt;
	}
} */

// dump hits
void dumpHits() {
	cout << "Dump hits" << endl;
	for (map<uint64_t, map<uint64_t, double>* >::iterator ref_it = hits.begin(), ref_eit = hits.end(); ref_it != ref_eit; ++ref_it) {
		for (map<uint64_t, double>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			cout << "Ref " << ref_it->first << " " << " Lease " << ri_it->first << " Hits " << (*hits[ref_it->first])[ri_it->first] << endl;
		}
		cout << endl;
	}
}

// dump costs
void dumpCosts() {
	cout << "Dump costs" << endl;
	for (map<uint64_t, map<uint64_t, double>* >::iterator ref_it = costs.begin(), ref_eit = costs.end(); ref_it != ref_eit; ++ref_it) {
		for (map<uint64_t, double>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			cout << "Ref " << ref_it->first << " " << " Lease " << ri_it->first << " Costs " << (*costs[ref_it->first])[ri_it->first] << endl;
		}
		cout << endl;
	}
}

// dump RI
void dumpRI() {
	uint64_t total_number_of_ri = 0;
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		std::set<uint64_t> riset;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			cout << "Ref " << ref_it->first << " RI " << ri_it->first << " CNT " << ri_it->second << endl;
			riset.insert(ri_it->first);
		}
		cout << "Ref " << ref_it->first << " RISETSIZE " << riset.size() << endl;
		total_number_of_ri += riset.size();
	}
	cout << "Average RISETSIZE for each reference " << double(total_number_of_ri) / RI.size() << endl;
}



void dumpLeases() {
	cout << "Dump leases" << endl;
	for (map<uint64_t, uint64_t>::iterator it = Lease.begin(), eit = Lease.end(); it != eit; ++it) {
		cout << it->first << " " << it->second << endl;
	}
	cout << endl;
}

// calculate hits and costs for each reuse interval 
void initHitsCosts() {

	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		hits[ref_it->first] = new map<uint64_t, double>;
		costs[ref_it->first] = new map<uint64_t, double>;
		(*hits[ref_it->first])[0] = 0;
		uint64_t total_hits = 0;
		(*costs[ref_it->first])[0] = 0;
		uint64_t total_cnt = 0;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			total_cnt += ri_it->second;
		}
		uint64_t pre_lease = 0;
		uint64_t pre_cost = 0;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			total_hits += ri_it->second;
			(*hits[ref_it->first])[ri_it->first] = total_hits;

			(*costs[ref_it->first])[ri_it->first] =  pre_cost + (ri_it->first - pre_lease) * total_cnt;
			total_cnt -= ri_it->second;
			pre_cost = (*costs[ref_it->first])[ri_it->first];
			pre_lease = ri_it->first;
		}
	}

}

// calculate PPUC
double getPPUC(uint64_t ref_id, uint64_t oldLease, uint64_t newLease) {
	//cout << "  getPPUC for ref " <<  ref_id << " oldLease " << oldLease << " newLease " << newLease << " ppuc: " << double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease])  << endl;

	//cout << (*hits[ref_id])[newLease] << " " << (*hits[ref_id])[oldLease] << " " << (*costs[ref_id])[newLease] << " " << (*costs[ref_id])[oldLease] << endl;
	//cout << (*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease] << " " << ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]) << endl;
	
	if (hits.find(ref_id) == hits.end() || costs.find(ref_id) == costs.end()) {
		cout << "No such ref for hits/costs" << endl;
		return -1;
	}
	if (hits[ref_id]->find(newLease) == hits[ref_id]->end() || costs[ref_id]->find(newLease) == costs[ref_id]->end()) {
		cout << "No RI/Newlease " << newLease << " for ref " << ref_id << endl;
		return -1;
	}
	
	if (hits[ref_id]->find(oldLease) == hits[ref_id]->end() || costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
		if (hits[ref_id]->find(oldLease) == hits[ref_id]->end()) {
			cout << "No hits for Oldlease " << oldLease << " for ref " << ref_id << endl;
		}
		if (costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
			cout << "No costs for Oldlease " << oldLease << " for ref " << ref_id << endl;
		}
		return -1;
	}

	return double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]);
}

// find max PPUC
void getMaxPPUC(bool*finished, uint64_t* ref_to_assign, uint64_t* newLease) {
	
	double maxPPUC = -1;
	uint64_t bestRef = -1;
	uint64_t bestLease = -1;

	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			//cout << "compare ri " << ri_it->first << " with lease " << Lease[ref_it->first] << " for reference " << ref_it->first << endl;
			if (ri_it->first > Lease[ref_it->first]) {
				double ppuc = getPPUC(ref_it->first, Lease[ref_it->first], ri_it->first);
				//cout << "  compare ppuc " << ppuc << " with max PPUC " << maxPPUC << endl;
				if (ppuc > maxPPUC) {
					maxPPUC = ppuc;
					bestRef = ref_it->first;
					bestLease = ri_it->first;
				}
				//break;
			}
		}
		//break;
	}

	if (maxPPUC != -1) {
		*finished = false;
		*ref_to_assign = bestRef;
		*newLease = bestLease;
	} else {
		*finished = true;
	}

	return;
}



// main OSL_ref alg
void OSL_ref(uint64_t CacheSize) {
	
	cout << "Start to init hits and costs" << endl;
	
	RIwithInfinite();
	initHitsCosts();

	//accessRatioCal();
	cout << "Finished to init hits and costs" << endl;
	
	//dumpHits();
	//dumpCosts();
	dumpRI();

	uint64_t N = refT;
	uint64_t totalCost = 0;
	uint64_t totalHits = 0;
	uint64_t targetCost = CacheSize * N;
	while(true) {
		bool finished = false;
		uint64_t ref_to_assign;
		uint64_t newLease;
		
		getMaxPPUC(&finished, &ref_to_assign, &newLease);
		//cout << "get max ppuc: " << finished << " " << ref_to_assign << " " << newLease << endl;

		if (finished == false) {
			totalCost += (*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]];
			totalHits += (*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]];
			Lease[ref_to_assign] = newLease;
			
			cout << "Assign lease " << newLease << " to ref " << ref_to_assign << " avg cache size " << double(totalCost) / N  << " miss ratio " << 1 - double(totalHits) / N << endl;
			
		} else {
			break;
		}
		
		if (totalCost < targetCost && targetCost != 0) {
			break;
		}
	}

	/*
	double totalCost = 0;
	double totalHitRatio = 0;
	double targetCost = CacheSize;
	while(true) {
		bool finished = false;
		uint64_t ref_to_assign;
		uint64_t newLease;
		
		getMaxPPUC(&finished, &ref_to_assign, &newLease);
		//cout << "get max ppuc: " << finished << " " << ref_to_assign << " " << newLease << endl;
		
		if (finished == false) {
			totalCost += ((*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
			totalHitRatio += ((*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
			Lease[ref_to_assign] = newLease;
			
			cout << "Assign lease " << newLease << " to ref " << ref_to_assign << " avg cache size " << totalCost << " miss ratio " << 1 - totalHitRatio << endl;
			
		} else {
			break;
		}
		
		if (totalCost < targetCost && targetCost != 0) {
			break;
		}
	}*/

	//dumpLeases();

	return;
}







