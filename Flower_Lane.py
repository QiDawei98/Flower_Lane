# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 20:15:20 2025

@author: qidaw
"""

import re
import numpy as np
from scipy.interpolate import RBFInterpolator
import math
import matplotlib.pyplot as plt
from multiprocessing import Pool
from multiprocessing import cpu_count
from numpy import linalg as LA
import os


def interpolator(args):
    xi, values, points = args  
    interp = RBFInterpolator(points, values, neighbors=60)
    return interp(xi)

if __name__ == '__main__':
    #use an integer smaller than cpu_count to resolve RAM issue
    cores = cpu_count()
    
    #Specify your path here
    file_path = r"G:\CHGCAR_mp-51"

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    scaling_factor = float(lines[1])

    m = re.search(r"^\s+ (-?\d+.\d+) \s+ (-?\d+.\d+) \s+ (-?\d+.\d+)", lines[2], re.VERBOSE) 
    a = np.array(m.groups(), dtype = float)
    a *= scaling_factor

    m = re.search(r"^\s+ (-?\d+.\d+) \s+ (-?\d+.\d+) \s+ (-?\d+.\d+)", lines[3], re.VERBOSE) 
    b = np.array(m.groups(), dtype = float)
    b *= scaling_factor

    m = re.search(r"^\s+ (-?\d+.\d+) \s+ (-?\d+.\d+) \s+ (-?\d+.\d+)", lines[4], re.VERBOSE) 
    c = np.array(m.groups(), dtype = float)
    c *= scaling_factor

    start_copying = False
    extracted_lines = []

    for line in lines[7:]:
        m = re.search(r"^\s* (\d\d+) \s+ (\d\d+) \s+ (\d\d+)$", line, re.VERBOSE)
        if m:
            NGXF = int(m.groups()[0])
            NGYF = int(m.groups()[1])
            NGZF = int(m.groups()[2])
            start_copying = True
            continue
        if re.search(r"\b augmentation \b", line, re.VERBOSE):
            break
        if start_copying:
            extracted_lines.append(line.strip())
        
    pristine_str = "\n".join(extracted_lines)
    pristine_data = np.array([float(val) for val in pristine_str.split()])

    n_points = NGXF * NGYF * NGZF
    assert pristine_data.shape == (n_points,)

    reshaped_data = pristine_data.reshape((NGXF, NGYF, NGZF), order='F')    

    x = []
    y = []
    z = []
    density = []

    for k in range(0, NGZF):    
        for j in range(0, NGYF):
            for i in range(0, NGXF):
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[i, j, k])
                
    assert len(x) == n_points
            
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in range(0, NGYF):
            for i in range(0, NGXF):
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])

    for k in range(0, NGZF):
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in range(NGXF):
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])

    for k in range(0, NGZF):
        for j in range(0, NGYF):
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])  
                
    for k in range(0, NGZF):
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF]) 
                
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in range(0, NGYF):
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])    
                    
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in range(0, NGXF):
                x.append(i /  NGXF * a[0] + j / NGYF * b[0] + k / NGZF * c[0])
                y.append(i /  NGXF * a[1] + j / NGYF * b[1] + k / NGZF * c[1])
                z.append(i /  NGXF * a[2] + j / NGYF * b[2] + k / NGZF * c[2])
                density.append(reshaped_data[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                
    
        
                     
            
    Vcell = abs(np.dot(a, np.cross(b, c)))   
    density = [value * 0.529177 ** 3 / Vcell for value in density] 

    hx = abs(a[0] / NGXF) + abs(b[0] / NGYF) + abs(c[0] / NGZF) * 0.8
    hy = abs(a[1] / NGXF) + abs(b[1] / NGYF) + abs(c[1] / NGZF) * 0.8
    hz = abs(a[2] / NGXF) + abs(b[2] / NGYF) + abs(c[2] / NGZF) * 0.8  

    xplus = [value + hx for value in x]
    xminus = [value - hx for value in x]
    yplus = [value + hy for value in y]
    yminus = [value - hy for value in y]
    zplus = [value + hz for value in z]
    zminus = [value - hz for value in z]

    points = np.column_stack((x, y, z))

    xi_combinations = [
        np.column_stack((xplus[:n_points], y[:n_points], z[:n_points])),
        np.column_stack((xminus[:n_points], y[:n_points], z[:n_points])),
        np.column_stack((x[:n_points], yplus[:n_points], z[:n_points])),
        np.column_stack((x[:n_points], yminus[:n_points], z[:n_points])),
        np.column_stack((x[:n_points], y[:n_points], zplus[:n_points])),
        np.column_stack((x[:n_points], y[:n_points], zminus[:n_points]))
    ]

    tasks = [(xi, density, points) for xi in xi_combinations]

    with Pool(cores if cores < 6 else 6) as pool:
        results = pool.map(interpolator, tasks)

    (density_xplus, density_xminus,
     density_yplus, density_yminus,
     density_zplus, density_zminus) = results
    
    fx = [(density_xplus[i] - density_xminus[i]) / (2 * hx) for i in range(n_points)]
    fy = [(density_yplus[i] - density_yminus[i]) / (2 * hy) for i in range(n_points)]
    fz = [(density_zplus[i] - density_zminus[i]) / (2 * hz) for i in range(n_points)]
        
    density_to_the_four_thirds = [abs(density[i]) ** (4 / 3) for i in range(n_points)]
    constant = 1 / (2 * (3 * np.pi**2)**(1/3))
    rdg = [constant * math.sqrt(fx[i]**2 + fy[i]**2 + fz[i]**2) / density_to_the_four_thirds[i] for i in range(n_points)]
        
    plt.scatter(density[:n_points], rdg)
    plt.ylim([0, 4])
    plt.xlim([0, 0.1])
    plt.show()
    
    reshaped_fx = np.array(fx).reshape((NGXF, NGYF, NGZF), order='F') 
    reshaped_fy = np.array(fy).reshape((NGXF, NGYF, NGZF), order='F')
    reshaped_fz = np.array(fz).reshape((NGXF, NGYF, NGZF), order='F')
    
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in range(0, NGYF):
            for i in range(0, NGXF):
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])

    for k in range(0, NGZF):
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in range(NGXF):
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])

    for k in range(0, NGZF):
        for j in range(0, NGYF):
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF]) 
                
    for k in range(0, NGZF):
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
    
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in range(0, NGYF):
            for i in [-1, -2, -3, -4, -5, NGXF, NGXF+1, NGXF+2, NGXF+3, NGXF+4]:
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF]) 
                    
    for k in [-1, -2, -3, -4, -5, NGZF, NGZF+1, NGZF+2, NGZF+3, NGZF+4]:
        for j in [-1, -2, -3, -4, -5, NGYF, NGYF+1, NGYF+2, NGYF+3, NGYF+4]:
            for i in range(0, NGXF):
                fx.append(reshaped_fx[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fy.append(reshaped_fy[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                fz.append(reshaped_fz[(i + NGXF) % NGXF, (j + NGYF) % NGYF, (k + NGZF) % NGZF])
                    

    
    tasks = []
    tasks.extend([(xi, fx, points) for xi in xi_combinations])
    tasks.extend([(xi, fy, points) for xi in xi_combinations])
    tasks.extend([(xi, fz, points) for xi in xi_combinations])
    
    with Pool(cores if cores < 18 else 18) as pool:
        results = pool.map(interpolator, tasks)
    
    (fx_xplus, fx_xminus,
     fx_yplus, fx_yminus,
     fx_zplus, fx_zminus,
     fy_xplus, fy_xminus,
     fy_yplus, fy_yminus,
     fy_zplus, fy_zminus,
     fz_xplus, fz_xminus,
     fz_yplus, fz_yminus,
     fz_zplus, fz_zminus) = results         
    
    fxx = [(fx_xplus[i] - fx_xminus[i]) / (2 * hx) for i in range(n_points)]
    fxy = [(fx_yplus[i] - fx_yminus[i]) / (2 * hy) for i in range(n_points)]
    fxz = [(fx_zplus[i] - fx_zminus[i]) / (2 * hz) for i in range(n_points)]

    fyx = [(fy_xplus[i] - fy_xminus[i]) / (2 * hx) for i in range(n_points)]
    fyy = [(fy_yplus[i] - fy_yminus[i]) / (2 * hy) for i in range(n_points)]
    fyz = [(fy_zplus[i] - fy_zminus[i]) / (2 * hz) for i in range(n_points)]
    
    fzx = [(fz_xplus[i] - fz_xminus[i]) / (2 * hx) for i in range(n_points)]
    fzy = [(fz_yplus[i] - fz_yminus[i]) / (2 * hy) for i in range(n_points)]
    fzz = [(fz_zplus[i] - fz_zminus[i]) / (2 * hz) for i in range(n_points)]
    
    Hessian = []
    for i in range(n_points):
        Hessian.append(np.array([[fxx[i], fyx[i], fzx[i]],
                                  [fxy[i], fyy[i], fzy[i]],
                                  [fxz[i], fyz[i], fzz[i]]]))
    for i in range(n_points):
       Hessian[i] = (Hessian[i] + Hessian[i].T) / 2 
    
    l2 = []
    count=0
    for i in range(n_points):
        eigenvalues = np.sort(LA.eigh(Hessian[i])[0])
        if not np.all(np.isreal(eigenvalues)):
            print(eigenvalues)
            count +=1
        assert np.all(np.isreal(eigenvalues))
        assert len(eigenvalues) == 3
        l2.append(eigenvalues[1])
        
    sl2rho = []
    for i in range(n_points):
        sl2rho.append(math.copysign(density[i], l2[i])) 
    
    rdg = np.array(rdg)
    sl2rho = np.array(sl2rho)
    density = np.array(density)
    rdg[density[:n_points] > 0.05] = 100

    plt.figure('NCI Scatter')
    index = np.where((rdg <= 2) & (abs(sl2rho) <= 0.05), True, False)
    x = sl2rho[index]
    y = rdg[index]
    plt.scatter(x, y,  c=x, cmap='jet', vmin = -0.04, vmax = 0.02)
    plt.colorbar()
    plt.xlabel(r'$sign(\lambda_2)\rho$ (a.u.)')
    plt.ylabel(r'$RDG$')
    plt.xlim([-0.05, 0.05])
    plt.ylim([0, 2])
    plt.show()
    
    for i in range(n_points):
        if rdg[i] > 100:
            rdg[i] = 100
    
    filename = os.path.basename(file_path)
    _, ext = os.path.splitext(filename)
    if not ext:
        filename += ".vasp"
    
    header_list = []
    for index, line in enumerate(lines):
        header_list.append(line)
        if index >= 7:
            if re.search(r"^\s*(\d\d+)\s+(\d\d+)\s+(\d\d+)\s*$", line, re.VERBOSE):
                break    
    header_content = "".join(header_list)
    
    rdgout = rdg * Vcell * 1.8897259886 ** 3

    with open("rdg_" + filename, 'w') as file:
        file.write(header_content)
        for i in range(0, len(rdgout), 5):
            line = '  '.join(f"{rdgout[j]:.11E}" for j in range(i, min(i + 5, len(rdgout))))
            file.write(line + '\n')
    
    filename = re.sub(r'chgcar', "ΨΗΓΨΑΡ", filename, flags=re.IGNORECASE | re.VERBOSE)
    for i in range(n_points):
        if sl2rho[i] < -0.04:
            sl2rho[i] = -0.04
        elif sl2rho[i] > 0.02:
            sl2rho[i] = 0.02
    with open("sl2rho_" + filename, 'w') as file:
        file.write(header_content)
        for i in range(0, len(sl2rho), 5):
            line = '  '.join(f"{sl2rho[j]:.11E}" for j in range(i, min(i + 5, len(sl2rho))))
            file.write(line + '\n')           