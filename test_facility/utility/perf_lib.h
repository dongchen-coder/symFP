#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

#define LLC_MISS_EVENT  0xd120
#define LLC_REF_EVENT   0xd104
#define MAX_PERF_EVENT_OPEN			3
#define MAX_PERF_READ_BUF_LENGTH	4096

using namespace std;

struct perf_event_attr pe;
int fds[MAX_PERF_EVENT_OPEN];
uint64_t ids[MAX_PERF_EVENT_OPEN];
char buf[MAX_PERF_READ_BUF_LENGTH];
struct read_format* rf = (struct read_format *) buf;
struct event_format* el[MAX_PERF_EVENT_OPEN];

struct read_format {
	uint64_t nr;
	struct {
		uint64_t value;
		uint64_t id;
	} values[];
};

struct event_format {
	string name;
	uint64_t raw_value;
	event_format(string name, uint64_t raw_value) {
		this->name = name;
		this->raw_value = raw_value;
	}
};

long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
               int cpu, int group_fd, unsigned long flags)
{
	int ret;

	ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
	            group_fd, flags);
	return ret;
}

int
read_event_from_file() {
	// Create a text string, which is used to output the text file
	string perf_event;
	// Read from the text file
	fstream pfevt_file;
	pfevt_file.open("PERF_EVENT_LIST.txt");
	int evt_cnt = 0;
	// Use a while loop together with the getline() function to read the file line by line
	string name;
	uint64_t raw_value;
	for(string word; pfevt_file >> word; ) {
		// Output the text from the file
		cout << word << endl;
		// this is a event raw value
		if (evt_cnt % 2 == 0) {
			name = word;
		} else {
			raw_value = stoi(word, 0, 16);
			el[evt_cnt / 2] = new event_format(name, raw_value);
		}
		if (evt_cnt / 2 >= MAX_PERF_EVENT_OPEN) {
			fprintf(stderr, "Support %d event at most. Only the first %d events will be traced\n", MAX_PERF_EVENT_OPEN, MAX_PERF_EVENT_OPEN);
			break;
		}
		evt_cnt++;
	}
	// Close the file
	pfevt_file.close();
	return (evt_cnt / 2);
}

void
perf_init() {
	int evt_cnt = read_event_from_file();
	// load event from file
	for (int i = 0; i < evt_cnt; i++) {
		memset(&pe, 0, sizeof(struct perf_event_attr));
		pe.type = PERF_TYPE_RAW;
		pe.size = sizeof(struct perf_event_attr);
		pe.config = el[i]->raw_value;
		pe.disabled = 1;
		pe.exclude_kernel = 1;
		pe.exclude_hv = 1;
		pe.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_ID;
		if (i == 0) {
			fds[i] = perf_event_open(&pe, 0, -1, -1, 0); 
		} else {
			fds[i] = perf_event_open(&pe, 0, -1, fds[i-1], 0); 
		}
		if (fds[i] == -1) {
			fprintf(stderr, "Error opening leader %llx\n", pe.config);
			exit(EXIT_FAILURE);
		}
		ioctl(fds[i], PERF_EVENT_IOC_ID, &ids[i]);
	}
}

void
perf_start() {
	ioctl(fds[0], PERF_EVENT_IOC_RESET, 0);
	ioctl(fds[0], PERF_EVENT_IOC_ENABLE, 0);
}

void 
perf_end() {
	ioctl(fds[0], PERF_EVENT_IOC_DISABLE, 0);
}


void
perf_print() {
	read(fds[0], &buf, sizeof(buf));
	for (int i = 0; i < rf->nr; i++) {
		for (int j = 0; j < MAX_PERF_EVENT_OPEN; j++) {
			if (rf->values[i].id == ids[j]) {
				printf("%s, %lld\n", el[i]->name.c_str(), rf->values[i].value);
				break;
			}
		}
	}
	for (int i = 0; i < MAX_PERF_EVENT_OPEN; i++) {
		close(fds[i]);
	}

}
