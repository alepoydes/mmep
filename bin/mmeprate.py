#!/usr/bin/env python3

import argparse
import numpy as np
import python.miser3 as miser

def translation_volumes(system, state):
	gen0=system.generator([1,0],state); v0=miser.norm(gen0)
	gen1=system.generator([0,1],state); v1=miser.norm(gen1)
	angle=np.arccos(np.sum(gen0*gen1)/v0/v1) if v0>0 and v1>0 else 0
	v01=v0*v1*np.sin(angle)
	return v0,v1,v01

class bcolors:
    FAIL = '' #'\033[91m'
    END = '' #'\033[0m'
    BOLD = '' #'\033[1m'
    UNDERLINE = '' #'\033[4m'

print('Transition rate calculator')
# parse command line
parser=argparse.ArgumentParser(
	description='Compute transition rates given MEP',
	epilog='Energy must be in meV, temperature in K.')
parser.add_argument('MEP', metavar='<mep.oct>', type=argparse.FileType('r'), help='Path to MEP saved by bin/mepx.')
parser.add_argument('T', metavar='{temperature}', type=float, nargs='+', help='Values of temperature')
parser.add_argument('--muB-factor', dest='muBfactor', default=1., type=float, help='Magnetic moment is assumed to be <muB-factor>*<Bohr magneton>')
parser.add_argument('-c', dest='color_output', action='store_true', help='Use color output')
args=parser.parse_args()

if args.color_output:
	bcolors.FAIL='\033[91m'
	bcolors.END = '\033[0m'
	bcolors.BOLD = '\033[1m'
	bcolors.UNDERLINE = '\033[4m'

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

print('Compiling energy function')
hessian=system.lambda_hessian()

print('Calculating eigenvalues for',bcolors.BOLD,'initial',bcolors.END,'state')
dataI=system.restricted_harmonic(path[0], hessian=hessian)
ergyI, peiI, neiI, cI, zeiI, zevI=dataI
if np.abs(ergyI-energy[0,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiI)
print('  Zero modes',zeiI)
print('  Positive modes',peiI[:3],'...')
v0,v1,v01=translation_volumes(system, path[0])
print('  Length of 1st translation mode',v0)
print('  Length of 2nd translation mode',v1)
print('  Area of translation plane',v01)


print('Calculating eigenvalues for',bcolors.BOLD,'transition',bcolors.END,'state')
dataS=system.restricted_harmonic(path[saddle], hessian=hessian)
ergyS, peiS, neiS, cS, zeiS, zevS=dataS
if np.abs(ergyS-energy[saddle,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiS)
print('  Zero modes',zeiS)
print('  Positive modes',peiS[:3],'...')
v0,v1,v01=translation_volumes(system, path[saddle])
print('  Length of 1st translation mode',v0)
print('  Length of 2nd translation mode',v1)
print('  Area of translation plane',v01)

print('Calculating eigenvalues for',bcolors.BOLD,'final',bcolors.END,'state')
dataE=system.restricted_harmonic(path[-1], hessian=hessian)
ergyE, peiE, neiE, cE, zeiE, zevE=dataE
if np.abs(ergyE-energy[-1,-1])>1e-8:
	print(bcolors.BOLD,'WARNING',bcolors.END,' stored energy is incorrect',sep='')
print('  Negative modes',neiE)
print('  Zero modes',zeiE)
print('  Positive modes',peiE[:3],'...')
v0,v1,v01=translation_volumes(system, path[-1])
print('  Length of 1st translation mode',v0)
print('  Length of 2nd translation mode',v1)
print('  Area of translation plane',v01)


print('Calculating rates without zero-modes and multiple saddles')
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