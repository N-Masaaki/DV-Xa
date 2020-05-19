#! /usr/bin/env python
# -*- encoding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import re

num = int(np.genfromtxt('l04', skip_header=1, max_rows=1)[0])
labelx_list = []
f = open('l04', 'r')
raw = f.readlines()
for n in range(1, len(raw)):
	c = raw[n].split()
	for i in range(len(c)):
		labelx_list.append(int(c[i]))
		if n > 1 and int(c[i]) == num:
			break
	if n > 1 and str(num) in c:
		break
labelx_list.pop(0)

f = open('d04', 'r')
raw = f.readlines()
for n in range(len(raw)):
	if re.match(r' \n', raw[n]):
		header = n
data_y = np.genfromtxt('d04', skip_header = header+1,delimiter='\t', dtype=str)
labely_list = []
for n in range(len(data_y)):
	for m in range(0, len(data_y[n].split()), 3):
		 labely_list.append(data_y[n].split()[m]+data_y[n].split()[m+1]+'_'+data_y[n].split()[m+2])

def lvlshm(A):
	plot_list = []
	if type(A) == list:
		for i in range(len(A)):
			for n in range(len(labely_list)):
				if re.match(r'%s[A-Z]+' %A[i], labely_list[n]):
					if labelx_list[n] > 0:
						plot_list.append(labelx_list[n])
	elif type(A) == str:			
		for n in range(len(labely_list)):
			if re.match(r'([0-9]+)%s' %A, labely_list[n]):
				if labelx_list[n] >0:
					plot_list.append(labelx_list[n])
	plot_list.append(num+1)
	
	data = pd.read_csv('l08', header=None)
	fig, axes = plt.subplots(1, len(plot_list), sharey=True)
	plot_index = 1
	for i in plot_list:
		m = i 
		step = (i - 1) // 40
		if 40 < i:
			i = i % 40
			if i == 0:
				i == 40
		x, y = [], []
		for n in range(0,len(data), math.ceil(num/40)):
			x.append(data[0][n])
			if step == 0:
				y.append(data[i][n])
			else:
				y.append(data[i][n+step])
		if m == num+1:
			x0, x1, y0, y1 = [], [], [], []
			for j in range(len(y)):
				if y[j] == 1:
					x1.append(x[j])
					y1.append(y[j])
				else:
					x0.append(x[j])
					y0.append(y[j] + 1)
			axes[0].barh(x1, y1, height=0.2, color='b')
			axes[0].barh(x0, y0, height=0.2, hatch='|||', color='white')
			axes[0].set_xticks(np.arange(0,1.4,0.2))
			axes[0].set_yticks(np.arange(-20,20,5))
			axes[0].set_ylim([-20,20])
			axes[0].spines['top'].set_visible(False)
			axes[0].spines['bottom'].set_visible(False)
			axes[0].spines['right'].set_visible(False)
			axes[0].set_xlabel('Total')
			axes[0].set_ylabel('Energy (eV)')
			axes[0].tick_params(bottom=False, labelbottom=False)
		else:
			axes[plot_index].barh(x, y, height=0.2, color='r')
			axes[plot_index].set_xticks(np.arange(0,1.4,0.2))
			axes[plot_index].set_yticks(np.arange(-20,20,5))
			axes[plot_index].set_ylim([-20,20])
			axes[plot_index].spines['top'].set_visible(False)
			axes[plot_index].spines['bottom'].set_visible(False)
			axes[plot_index].spines['right'].set_visible(False)
			axes[plot_index].spines['left'].set_visible(False)
			axes[plot_index].set_xlabel(labely_list[labelx_list.index(m)])
			axes[plot_index].tick_params(bottom=False, left=False, labelbottom=False)
			plot_index += 1
	fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0,hspace=0)
	fig.tight_layout()
	fig.show()
	return plot_list
			
def dos_check(A):
	plot_list = []
	if type(A) == list:
		for i in range(len(A)):
			for n in range(len(labely_list)):
				if re.match(r'%s[A-Z]+' %A[i], labely_list[n]):
					if labelx_list[n] > 0:
						plot_list.append(labelx_list[n])
	elif type(A) == str:			
		for n in range(len(labely_list)):
			if re.match(r'([0-9]+)%s' %A, labely_list[n]):
				if labelx_list[n] >0:
					plot_list.append(labelx_list[n])
	plot_list.append(num+1)
	
	data = pd.read_csv('dd7', header=None)
	plt.figure()
	for i in plot_list:
	#for i in range(1,num+2):
		m = i
		step = (i - 1) // 30
		if 30 < i:
			i = i % 30
			if i == 0:
				i = 30
		x, y = [], []
		for n in range(0,len(data), math.ceil(num/30)):
			x.append(data[0][n])
			if step == 0:
				y.append(data[i][n])
			else:
				y.append(data[i][n+step])
		if m <= num:
			plt.plot(y,x, label=labely_list[labelx_list.index(m)])
			#plt.plot(y,x)
		else:
			plt.plot(y,x, label='all', color='k')
	plt.ylabel('Energy (eV)')
	#plt.yticks(np.arange(-3,7,1))
	#plt.ylim([-3,6])
	plt.legend(prop={'size': 8}, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, ncol=2)
	plt.tight_layout()
	plt.show()


def dos_check2(A):
	plot_list = []
	if type(A) == list:
		for i in range(len(A)):
			for n in range(len(labely_list)):
				if re.match(r'%s[A-Z]+' %A[i], labely_list[n]):
					if labelx_list[n] > 0:
						plot_list.append(labelx_list[n])
	elif type(A) == str:			
		for n in range(len(labely_list)):
			if re.match(r'([0-9]+)%s' %A, labely_list[n]):
				if labelx_list[n] >0:
					plot_list.append(labelx_list[n])
	#plot_list.append(num+1)
	
	data = pd.read_csv('dd7', header=None)
	plt.figure(figsize=(3,4))
	for i in plot_list:
	#for i in range(1,num+2):
		m = i
		step = (i - 1) // 30
		if 30 < i:
			i = i % 30
			if i == 0:
				i = 30
		x, y = [], []
		for n in range(0,len(data), math.ceil(num/30)):
			x.append(data[0][n])
			if step == 0:
				y.append(data[i][n])
			else:
				y.append(data[i][n+step])
		if m <= num:
			plt.plot(y,x)#, label=labely_list[labelx_list.index(m)])
			#plt.plot(y,x)
		else:
			plt.plot(y,x, label='all')
	plt.ylabel('Energy (eV)')
	plt.xticks(np.arange(0,2.51,0.5))
	plt.xlim([0,2.51])
	#plt.yticks(np.arange(-3,7,1))
	#plt.ylim([-3,6])
	#plt.legend(prop={'size': 8}, loc='upper right', bbox_to_anchor=(1, 1), borderaxespad=0, ncol=2)
	plt.tight_layout()
	plt.show()









