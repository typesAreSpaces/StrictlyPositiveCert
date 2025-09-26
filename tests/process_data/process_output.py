#!/usr/bin/env python3

# This programs parses and write rows
# of the following format
# ```
# Number of isolated points & Max degree & Average degree & Average time &
# Total time & Ratio of success
# ```

import re
import sys

NUM_DIGITS_ROUND = 3

STATE_READ_NORMAL = 0
STATE_READ_BASIS  = 1
STATE_READ_TIME   = 2
STATE_READ_DEGREE = 3

def average(l):
    sum = 0
    size = len(l)
    i = 0
    while True:
        if i >= size:
            break
        sum += l[i]
        i = i + 1
    return round(float(sum)/size, NUM_DIGITS_ROUND)

def _median(l):
    if (len(l) == 0):
        return "-"
    _l = sorted(l)
    return _l[int(len(l)/2)]

class State:
    def __init__(self):
        self.reset()

    def readBasis(self, line):
        self.num_isolated_points = len(line.split(',')[0].split('*'))
        self.num_tests += 1

    def readTime(self, line):
        self.times.append(float(line))
        self.num_solved += 1

    def readDegree(self, line):
        self.degrees.append(int(line))

    def reset(self):
        self.times               = []
        self.degrees             = []
        self.num_isolated_points = 0
        self.num_solved          = 0
        self.num_tests           = 0

    def print1(self):
        if len(self.times) > 0:
            print(
                    self.num_isolated_points, "&", 
                    max(self.degrees), "&", 
                    average(self.degrees), "&",
                    _median(self.degrees), "&",
                    average(self.times), "&", 
                    _median(self.times), "&",
                    round(sum(self.times), NUM_DIGITS_ROUND), "&", 
                    str(self.num_solved) + '/' + str(self.num_tests), "\\\\") 
        elif self.num_tests > 0:
            print(
                    "-", "&", 
                    "-", "&", 
                    "-", "&",
                    "-", "&",
                    "-", "&", 
                    "-", "&", 
                    "-", "\\\\") 

    def print2(self):
        if len(self.times) > 0:
            print(
                    self.degrees[0], "&", 
                    self.times[0], "\\\\") 
        elif self.num_tests > 0:
            print(
                    "-", "&", 
                    "-", "\\\\") 

class Entries:
    def __init__(self, num_entries):
        self.state = STATE_READ_NORMAL
        self.curr_entry = -1
        self.num_entries = num_entries
        self.entries = [State() for x in range(num_entries)]

    def next(self):
        self.curr_entry = (self.curr_entry + 1) % self.num_entries

    def readBasis(self, line):
        self.state = STATE_READ_NORMAL
        self.entries[self.curr_entry].readBasis(line)

    def readTime(self, line):
        self.state = STATE_READ_NORMAL
        self.entries[self.curr_entry].readTime(line)

    def readDegree(self, line):
        self.state = STATE_READ_NORMAL
        self.entries[self.curr_entry].readDegree(line)

    def reset(self):
        self.state = STATE_READ_NORMAL
        [x.reset() for x in self.entries]

    def print1(self):
        [x.print1() for x in self.entries]

    def print2(self):
        [x.print2() for x in self.entries]

if __name__ == "__main__":
    line_matcher = re.compile('.*>>.*')
    file = open(sys.argv[1], 'r')
    num_entries = int(sys.argv[2])
    print_method = int(sys.argv[3])

    entries = Entries(num_entries)

    # Skip until ">> Start Benchmarks"
    while True:
        line = file.readline().strip()
        if line == ">> Start Benchmarks":
            break

    while True:
        line = file.readline()

        if not line:
            break
 
        if entries.state == STATE_READ_NORMAL and line_matcher.match(line):
            line = line.strip()
            if line == ">> Start benchmark":
                # Write row from previous benchmark results
                if print_method == 0:
                    entries.print1()
                else:
                    entries.print2()
                # Reset state
                entries.reset() 
            if line == ">> Test":
                entries.next()
            if line == ">> basis":
                entries.state = STATE_READ_BASIS
            if line == ">> Time taken":
                entries.state = STATE_READ_TIME
            if line == ">> Degree size":
                entries.state = STATE_READ_DEGREE
            continue

        if entries.state == STATE_READ_BASIS:
            entries.readBasis(line)
            continue

        if entries.state == STATE_READ_TIME:
            entries.readTime(line)
            continue

        if entries.state == STATE_READ_DEGREE:
            entries.readDegree(line)
            continue

    # Write last row
    if print_method == 0:
        entries.print1()
    else:
        entries.print2()
