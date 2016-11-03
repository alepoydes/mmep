#!/usr/bin/env python3

import argparse
import numpy as np
import python.miser3 as miser

class bcolors:
    FAIL = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print('Transition rate calculator')
# parse command line
parser=argparse.ArgumentParser(
	description='Compute transition rates given MEP',
	epilog='Energy must be in meV, temperature in K.')
parser.add_argument('MEP', metavar='<mep.oct>', type=argparse.FileType('r'), help='Path to MEP saved by bin/mepx.')
parser.add_argument('T', metavar='{temperature}', type=float, nargs='+', help='Values of temperature')
parser.add_argument('--muB-factor', dest='muBfactor', default=1., type=float, help='Magnetic moment is assumed to be <muB-factor>*<Bohr magneton>')
args=parser.parse_args()

print(bcolors.BOLD,'Temperatures ',bcolors.END,args.T,sep='')
# Магнетон Бора
muB=5.7883818066E-2 # meV/T
# Постоянная Планка
hbar=6.582119514E-13 # meV*s
# Магнитный момент
mu=args.muBfactor*muB
# Множитель Ланде
ge = 2.00231930436
# Гиромагнитное отношение
gamma=ge*muB/hbar
# Постоянный множитель для предэкспоненты
gammaovermu=gamma/mu
# Константа Больцмана
kB=8.6173303E-2 # meV/K

print('Reading file')
system, path, energy, distance=miser.load_results(args.MEP)

saddle=np.argmax(energy[:,-1]) # index of saddle point
print(bcolors.BOLD,'Forward barrier',bcolors.END,' {:e}'.format(energy[saddle,-1]-energy[0,-1]),sep='')
print(bcolors.BOLD,'Backward barrier',bcolors.END,' {:e}'.format(energy[saddle,-1]-energy[-1,-1]),sep='')

print('Compilling energy function')
hessian=system.lambda_hessian()

print('Calculating eigenvalues for initial state')
dataI=system.restricted_harmonic(path[0], hessian=hessian)
ergyI, peiI, neiI, cI, zeiI, zevI=dataI
if np.abs(ergyI-energy[0,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiI)
print('  Zero modes',zeiI)
print('  Positive modes',peiI[:3],'...')

print('Calculating eigenvalues for transition state')
dataS=system.restricted_harmonic(path[saddle], hessian=hessian)
ergyS, peiS, neiS, cS, zeiS, zevS=dataS
if np.abs(ergyS-energy[saddle,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiS)
print('  Zero modes',zeiS)
print('  Positive modes',peiS[:3],'...')

print('Calculating eigenvalues for final state')
dataE=system.restricted_harmonic(path[-1], hessian=hessian)
ergyE, peiE, neiE, cE, zeiE, zevE=dataE
if np.abs(ergyE-energy[-1,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiE)
print('  Zero modes',zeiE)
print('  Positive modes',peiE[:3],'...')

print('Calculating rates')
rateI=[]; rateE=[];
for T in args.T:
	rateI.append(miser.rate(dataI, dataS, gammaovermu=gammaovermu, kT=T*kB))
	rateE.append(miser.rate(dataE, dataS, gammaovermu=gammaovermu, kT=T*kB))
rateI=np.array(rateI)	
rateE=np.array(rateE)

print(bcolors.BOLD,'Forward pre-exponents ',bcolors.END,np.prod(rateI[:,:-1],axis=1),sep='')
print(bcolors.BOLD,'Forward exponents ',bcolors.END,rateI[:,-1],sep='')
print(bcolors.BOLD,'Forward transition rates ',bcolors.END,np.prod(rateI,axis=1),sep='')

print(bcolors.BOLD,'Backward pre-exponents ',bcolors.END,np.prod(rateE[:,:-1],axis=1),sep='')
print(bcolors.BOLD,'Backward exponents ',bcolors.END,rateE[:,-1],sep='')
print(bcolors.BOLD,'Backward transition rates ',bcolors.END,np.prod(rateE,axis=1),sep='')