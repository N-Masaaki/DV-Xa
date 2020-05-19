#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import re
import periodic_table

cfg = "B10O21Pb6.xyz"
#cfg = "mp-1192421_B2Pb2O5.xyz"
xyz = (np.genfromtxt(cfg,skip_header=2,usecols=(1,2,3)))
#xyz = (np.genfromtxt(cfg,skip_header=2,usecols=(1,2,3)))
bohr_au = 0.529177      #ボーア半径原子単位(atomic unit) Å

Pb = xyz[0:384]
B = xyz[384:1024]
O = xyz[1024:]
#Pb = xyz[0:405]
#B = xyz[405:789]
#O = xyz[789:]
type = 3
         
def makefile(A, C_size, M_size):
    ex_xyz = xyz
    ex_Pb = Pb
    ex_B = B
    ex_O = O
    data_xyz = ex_xyz - ex_xyz[A]
    data_Pb = ex_Pb - ex_xyz[A]
    data_B = ex_B - ex_xyz[A]
    data_O = ex_O - ex_xyz[A]           #全座標を中心原子分移動
    f_c = open('f01', 'a+')
    f_m = open('f03', 'a+')
    f_c.write("| Z ||NEQ||      XD     ||      YD     ||      ZD     |\n")
    f_m.write('   2.00000   3.00000  -2.00000\n        1.0000000000        1.0000000000        1.0000000000\n')
    NEQ = 1
    NEQ_list, CHG_list = [], []
    MP = 0
    for n in range (len(data_xyz)):
        text_c = "{:>5}{:>5}{:>15.5f}{:>15.5f}{:>15.5f}\n"
        text_m = "{:>20.10f}{:>20.10f}{:>20.10f}{:>5}\n"
        if n == A:
            NEQ_list.append(NEQ)
            if data_xyz[n] in data_Pb:
                CHG_list.append(2.00000)
                f_c.write(text_c.format(82, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 1))
                NEQ += 1
                MP += 1
            elif data_xyz[n] in data_B: 
                CHG_list.append(3.00000)
                f_c.write(text_c.format(5, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 2))
                NEQ += 1
                MP += 1
            elif data_xyz[n] in data_O: 
                CHG_list.append(-2.00000)
                f_c.write(text_c.format(8, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 3))
                NEQ += 1
                MP += 1
        else:
            r = np.linalg.norm(data_xyz[n])
            if r <= C_size:
                NEQ_list.append(NEQ)
                if data_xyz[n] in data_Pb:
                    CHG_list.append(2.00000)
                    f_c.write(text_c.format(82, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                    NEQ += 1
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 1))
                    MP += 1
                elif data_xyz[n] in data_B:
                    CHG_list.append(3.00000) 
                    f_c.write(text_c.format(5, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                    NEQ += 1
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 2))
                    MP += 1
                elif data_xyz[n] in data_O:
                    CHG_list.append(-2.00000)
                    f_c.write(text_c.format(8, NEQ, data_xyz[n][0], data_xyz[n][1], data_xyz[n][2]))
                    NEQ += 1
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 3))
                    MP += 1
            elif C_size < r <= M_size:
                if data_xyz[n] in data_Pb:
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 1))
                    MP += 1
                elif data_xyz[n] in data_B:
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 2))
                    MP += 1
                elif data_xyz[n] in data_O:
                    f_m.write(text_m.format(data_xyz[n][0]/bohr_au, data_xyz[n][1]/bohr_au, data_xyz[n][2]/bohr_au, 3))
                    MP += 1
        print(n)
    f_c.write("---------------------------------------------\n|NEQ||  CHG   ||U/D||   RD   ||   VD   |    1\n")
    for i in range(len(NEQ_list)):
        text = "{:>5}{:>10.5f}\n"
        f_c.write(text.format(NEQ_list[i], CHG_list[i]))  
    f_c.write('''---------------------------------------------
    0     Unit     (0:angstrom  1:atomic)
    0     Spin     (0:non-spin  1:spin  )
    1     M.P.     (0:No        1:Yes   )
    0     Sample Point ( <100000, =0 autoset )''' )
    f_m.write('''    0    0    0
             0.00000             0.00000             0.00000
          9999.99999          9999.99999          9999.99999''')    
    f_c.close()
    f_m.close()
    text = "  {}    {}    0\n"
    with open('f03') as f_m:
        data = f_m.readlines()
    data.insert(0, text.format(MP, type))
    with open('f03', mode='w') as f_m:
        f_m.writelines(data)
        
    return NEQ_list, CHG_list
    
def DVfile_to_xyzfile():
    f01 = np.genfromtxt('f01', delimiter='\n', dtype=str)
    for i in range(len(f01)):
        if re.match(r'-+', f01[i]):
            f01_row = i - 1
            break
    f01_element = np.genfromtxt('f01', skip_header=1, max_rows=f01_row, usecols=(0), dtype=int)
    f01_xyz = np.genfromtxt('f01', skip_header=1, max_rows=f01_row, usecols=(2,3,4), dtype=float)
    fc = open('f01.xyz','a+')
    header_text = '{}\nPb{}B{}O{}\n'
    xyz_text = '{}\t{:.5f}\t{:.5f}\t{:.5f}\n'
    fc.write(header_text.format(f01_row, np.sum(f01_element==82), np.sum(f01_element==5), np.sum(f01_element==8)))
    for n in range(len(f01_xyz)):
        fc.write(xyz_text.format(periodic_table.element_dict(int(f01_element[n])), f01_xyz[n][0], f01_xyz[n][1], f01_xyz[n][2]))
    fc.close()
    
    f03_row = int(np.genfromtxt('f03',max_rows=1, usecols=(0)))
    f03_element = np.genfromtxt('f03', skip_header=3, max_rows=f03_row, usecols=(3), dtype=int)
    f03_xyz = np.genfromtxt('f03', skip_header=3, max_rows=f03_row, usecols=(0,1,2), dtype=float)
    adjust = [82,5,8]
    fm = open('f03.xyz','a+')
    header_text = '{}\nPb{}B{}O{}\n'
    xyz_text = '{}\t{:.10f}\t{:.10f}\t{:.10f}\n'
    fm.write(header_text.format(f03_row, np.sum(f03_element==1), np.sum(f03_element==2), np.sum(f03_element==3)))
    for n in range(len(f03_xyz)):
        fm.write(xyz_text.format(periodic_table.element_dict(adjust[f03_element[n]-1]), f03_xyz[n][0] * bohr_au, f03_xyz[n][1] * bohr_au, f03_xyz[n][2] * bohr_au))
    fm.close()






    
    

