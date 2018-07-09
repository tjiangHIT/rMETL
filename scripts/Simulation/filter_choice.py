#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  call_TE.py
 * @Package: argparse, sys, os, logging, pysam, Bio
 * @Description: Establish the ME callset
 * @author: tjiang
 * @date: Apr 24 2018
 * @version V1.0     
'''
import sys

def acquire_count_max(_list_):
	c = Counter(_list_)
	return c.most_common(1)[0][0]

flag_dic = {0:1, 16:2, 2048:0, 2064:0, 4:0}

STRAND = {'1':'+', '2':'-'}

cluster_dic = {}

strand_dic = {1:'+', 2:'-'}

def parse_name_tp(line):
	seq = line.split('*')
	return seq[4], seq[0]

def call(path):
	id_dic = dict()
	AlignmentFile = open(path, 'r')
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		local_info = parse_name_tp(seq[0])

		sub_type = seq[2]
		MAPQ = int(seq[4])
		if sub_type != "*" and MAPQ >= 20:
			if local_info[0] not in id_dic and local_info[1] == 'chr1':
				id_dic[local_info[0]] = 0
	AlignmentFile.close()
	for k in id_dic:
		print k

def run():
    call(sys.argv[1])

if __name__ == '__main__':
	run()
