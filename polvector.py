#!/usr/bin/python

#=======================================================================#
#  																# 
#  Oliver Gindele, UCL, 2014													#
#  Version 2.7													    		#		
# 																#									
#  read DLPOLY HISTORY, calculate local dipoles, total polarization and generate gnuplot output			 	#
#  NEW: correlation functions				 									#
#  use polmovie to create gif movies											    	#
#  built for perovsiktes												    		#
#=======================================================================#

import numpy as np
import sys,string, os, math
import argparse
import multiprocessing
from multiprocessing import Pool as Pool

version='2.7' 

#np.set_printoptions(threshold=np.nan)		#print whole array
np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=150)

#get command arguments
parser = argparse.ArgumentParser(description='Calculate local dipole moments from DLPOLY HISTORY file.')
parser.add_argument("dim",  type=int, nargs=3,  help="You must specify system dimensions (a b c)")
parser.add_argument('--file', '-f', type=str,  help="File name (default=HISTORY)", default='HISTORY', dest="fileName") 
parser.add_argument('--timestep', '-t', type=int,  help="select specific frame (timestep) to read", default=1, dest="frame")
parser.add_argument('--movie', '-m', type=int, nargs=3, help="select start timestep, end timestep and step to read multiple timestep in order to create a movie", dest="movie")
parser.add_argument('--average', '-a', type=str,  help="select direction (X,Y or Z) along polarization averages should be calculated (default = N, use N for no averaging)", 
default='N', dest="averageDirection")
parser.add_argument('--second average', '-s', type=str,  help="select second direction along polarization averages should be calculated", 
dest="secondAverageDirection")
parser.add_argument('--timeAverage',  help="Average over all selected frames (needed for calculation of piezoelectric coefficients)", action="store_true")
parser.add_argument('--densvar', '-d' ,type=float, help="Specify  a densvar factor that increases the amount of atoms each cell holds (default=1). Increase if non defined polarization", 
default=1.0, dest="densvar")
parser.add_argument('--binFactor', '-b' ,type=int, nargs=3, help=("Specify an integer factor to increase the cell size for the polarization binning. This will average polarization vectors over multiple cells for better visualization (default=1 1 1)"), 
default=[1 ,1 ,1], dest="binFactor")
parser.add_argument('--cation', '-c' ,type=str, help=("Select, if only one species (Ti or Zr) should be respected to calculate polarization (default = Ti,Zr)"), 
default='default', dest="cationSpecify")
parser.add_argument('--coreNumber', '-n' ,type=int, help=("Number of cores"), 
default='1', dest="coreNumber")
parser.add_argument('--correlation', help="Calculate spatial correlation of observables)", action="store_true")
parser.add_argument('--distribution', type=float, nargs=3,  help="Calculate polarization distribution. Print out Positions and neighbour distribution of ions in selected polarization range (e.g. 30, 50) -> use 1to switch on (third integer)")

args = parser.parse_args()


###########################################################
def getInfo(inputFile):
	"get info about MD run from HISTORY file"	
										
	inputFile=open(inputFile,'r')

	title=inputFile.readline()                                          	#read in header

	print "***** P O L V E C T O R *****" ,'\n'
	print "Version: " ,version
	print "Number of cores: ",args.coreNumber
	print "\n", "Prosessing System:", title

	line=inputFile.readline()                                          	#read in first line
	totalAtom=int(string.split(line)[2]+'\n')             	#get number of atoms 
	totalTimestep=int(string.split(line)[3]+'\n')	 	#get number of at timestep
	return totalAtom,totalTimestep

def getCell(inputFile,frame):
	"get Cell Vectors  HISTORY file"
	inputFile=open(inputFile,'r')
        line=inputFile.readline()  
        
        while(line!=""): 
        	if string.split(line)[0][:8]=='timestep' :
			line=inputFile.readline()
			cellVectorA=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])
			line=inputFile.readline()  
			cellVectorB=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])
			line=inputFile.readline()  
			cellVectorC=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])		
			break
		line=inputFile.readline()  						
	return cellVectorA, cellVectorB, cellVectorC
	
def getPos(inputFile,cellVectorA, cellVectorB, cellVectorC,totalAtom,frame,superCell):   
	"get positions from HISTORY file and save in array"
	atomList=['Pb', 'Ti', 'Zr','Mg','Nb', 'O','Pbs', 'Tis', 'Zrs','Os','Mgs','Nbs']
        position=np.empty(shape=(totalAtom,4),dtype=object)

        i=0                                                                    		 					#counter for position
        inputFile=open(inputFile,'r')
        line=inputFile.readline()  
	
        while(line!=""): 
       		line=inputFile.readline()
        	while(line!="") and i<totalAtom:
			if string.split(line)[0][:3] in atomList:                       
		        	position[i,0]=str(string.split(line)[0][:3])             				 #get atom name
		                        
		                line=inputFile.readline()                                      				 #advance a line
		                      
		                position[i,1]=float(string.split(line)[0])%cellVectorA[0]                  #get x pos
		                position[i,2]=float(string.split(line)[1])%cellVectorB[1]                  #get y pos
		                position[i,3]=float(string.split(line)[2])%cellVectorC[2]                  #get z pos
		                i+=1
			line=inputFile.readline()
		line=inputFile.readline() 
	
        inputFile.close()
        return position

def getAverageCellVectors(inputFile,cellVectorA, cellVectorB, cellVectorC):
	"get average cell vectors from OUTPUT to calculate strain accurately"
	foundLine=0; outputTrue=0
	
	#check if OUTPUT exists. If yes, read averaged cell vectors from there.
	if os.path.isfile(inputFile):
		searchFile=open(inputFile,'r')
		for line in searchFile:
			if line.startswith('                Average cell vectors'):
				line=searchFile.next();line=searchFile.next()		#skip lines
				averageCellVectorA=float(line.split()[0]); line=searchFile.next()
				averageCellVectorB=float(line.split()[1]); line=searchFile.next()
				averageCellVectorC=float(line.split()[2]); 
				foundLine=1; outputTrue=1
				break
				
	#if no OUTPUT, use the cell parameters at the time step (only ok for local strain)
	else:
		averageCellVectorA=cellVectorA[0]; averageCellVectorB=cellVectorB[1]; averageCellVectorC=cellVectorC[2];
	#if OUTPUT exists but run not completed
	if os.path.isfile(inputFile) and foundLine==0: 
		averageCellVectorA=cellVectorA[0]; averageCellVectorB=cellVectorB[1]; averageCellVectorC=cellVectorC[2];

	return averageCellVectorA, averageCellVectorB, averageCellVectorC, outputTrue
	
       
def increasePos(posPeriod, cellVectorA, cellVectorB, cellVectorC, binCounter, superCell):
	"increase position array and bincounter array by two layers in each direction"
	posPeriod=np.copy(posPeriod)
	posPeriodLarge=np.empty(shape=(superCell[0]+2,superCell[1]+2,superCell[2]+2,20*densvar,4),dtype=object)
	binCounterLarge=np.zeros(shape=(superCell[0]+2,superCell[1]+2,superCell[2]+2),dtype=int)

	#print "Increasing system size"

	posPeriodLarge[1:superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1]=posPeriod 					#copy position array
	binCounterLarge[1:superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1]=binCounter 					#copy binCount array

	#along X
	posPeriodLarge[0,1:superCell[1]+1,1:superCell[2]+1]=posPeriod[superCell[0]-1,:,:];  posPeriodLarge[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1]=posPeriod[0,:,:]; 
	binCounterLarge[0,1:superCell[1]+1,1:superCell[2]+1]=binCounter[superCell[0]-1,:,:];  binCounterLarge[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1]=binCounter[0,:,:];
	
	for j in xrange(superCell[1]+2):										#shift atom postions by one cell vector accordingly
			for k in xrange(superCell[2]+2): 	
				for ion in xrange(binCounterLarge[0,j,k]):
					posPeriodLarge[0,j,k,ion,1]-=cellVectorA[0]	
				for ion in xrange(binCounterLarge[superCell[0]+1,j,k]):
					posPeriodLarge[superCell[0]+1,j,k,ion,1]+=cellVectorA[0]	
	#along Y
	posPeriodLarge[:,0,1:superCell[2]+1]=posPeriodLarge[:,superCell[1],1:superCell[2]+1];  posPeriodLarge[:,superCell[1]+1,1:superCell[2]+1]=posPeriodLarge[:,1,1:superCell[2]+1]; 
	binCounterLarge[:,0,1:superCell[2]+1]=binCounterLarge[:,superCell[1],1:superCell[2]+1];  binCounterLarge[:,superCell[1]+1,1:superCell[2]+1]=binCounterLarge[:,1,1:superCell[2]+1]; 

	for i in xrange(superCell[0]+2):
			for k in xrange(superCell[2]+2): 	
				for ion in xrange(binCounterLarge[i,0,k]):
					posPeriodLarge[i,0,k,ion,2]-=cellVectorB[1]	
				for ion in xrange(binCounterLarge[i,superCell[1]+1,k]):
					posPeriodLarge[i,superCell[1]+1,k,ion,2]+=cellVectorB[1]
	#along Z
	posPeriodLarge[:,:,0]=posPeriodLarge[:,:,superCell[2]];  posPeriodLarge[:,:,superCell[2]+1]=posPeriodLarge[:,:,1]; 
	binCounterLarge[:,:,0]=binCounterLarge[:,:,superCell[2]];  binCounterLarge[:,:,superCell[2]+1]=binCounterLarge[:,:,1]; 

	for i in xrange(superCell[0]+2):
			for j in xrange(superCell[1]+2): 	
				for ion in xrange(binCounterLarge[i,j,0]):
					posPeriodLarge[i,j,0,ion,3]-=cellVectorC[2]	
				for ion in xrange(binCounterLarge[i,j,superCell[2]+1]):
					posPeriodLarge[i,j,superCell[2]+1,ion,3]+=cellVectorC[2]
	
	return posPeriodLarge,binCounterLarge

def calcPol(position,cellVectorA, cellVectorB, cellVectorC,frame,superCell,cationSpecify,nCell):        
	"calculate local polarization centered on B cation"
	polarization=np.zeros(shape=(nCell,6)); PbList=np.empty([8,3]); OList=np.empty([6,3])
	strain=np.zeros(shape=(nCell,9)); cellStrain=np.zeros(shape=(1,6))
	tetFactor=np.zeros(shape=(nCell,4)); unitVolume=np.zeros(shape=(nCell,7))
	octVector=np.zeros(shape=(nCell,6)); octAngle=np.zeros(shape=(nCell,6))
	cationPosition=np.empty(shape=(nCell,4))
	
	#ionic charges from force-field (PZT, O. Gindele 2015)
	PbCharge=11.697;	 	PbsCharge=-10.224
	TiCharge= 7.997;  	 	TisCharge= -5.047
	ZrChrage=11.229;   	ZrsCharge=-8.279
	OCharge= 2.467;     	OsCharge =-3.941
	cellVolume=cellVectorA[0]*cellVectorB[1]*cellVectorC[2]/(superCell[0]* superCell[1]* superCell[2])  #calculate cell volume

 	#bin values in order to make system a 3d array
	sizeA=np.linalg.norm(cellVectorA)/superCell[0]; sizeB=np.linalg.norm(cellVectorB)/superCell[1]; sizeC=np.linalg.norm(cellVectorC)/superCell[2]
	posPeriod=np.empty(shape=(superCell[0],superCell[1],superCell[2],20*densvar,4),dtype=object)
	binCounter=np.zeros(shape=(superCell[0],superCell[1],superCell[2]),dtype=int)				#Counter for number of atoms in each bin

	#print "Binning position values"	
	for m in xrange(len(position)):														#loop over all polarization values
															
		binX=math.floor(position[m,1]/sizeA); binY=math.floor(position[m,2]/sizeB); binZ=math.floor(position[m,3]/sizeC)

		posPeriod[binX,binY,binZ,binCounter[binX,binY,binZ],0]=position[m,0]
		posPeriod[binX,binY,binZ,binCounter[binX,binY,binZ],1:4]=position[m,1:4]
		binCounter[binX,binY,binZ]+=1
		
	#increase system
	posPeriod,binCounter=increasePos(posPeriod,cellVectorA, cellVectorB, cellVectorC,binCounter,superCell)

	#if only one cation selected change list accordingly:
	if cationSpecify=='Ti' or cationSpecify=='ti':
		cationList=['Ti']
	elif cationSpecify=='Zr' or cationSpecify=='zr':
		cationList=['Zr']
	elif cationSpecify=='default':
		cationList=['Ti','Zr']

	#get average cell vectors to calculate piezoelectric tensor:
	averageCellVectorA,averageCellVectorB,averageCellVectorC,outputTrue=getAverageCellVectors('OUTPUT',cellVectorA, cellVectorB, cellVectorC)
	
	m=0													#set counter for B cations	
	#print "Calculating polarization"
	for i in xrange(1,superCell[0]+1):
		for j in xrange(1,superCell[1]+1):
			for k in xrange(1,superCell[2]+1): 			
				
				for ii in xrange(binCounter[i,j,k]):
					cationPosition[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4]); 
					
					if str(posPeriod[i,j,k,ii,0])=='Ti':			#set array element 0 for Ti
						cationPosition[m,3]=0
					elif str(posPeriod[i,j,k,ii,0])=='Zr':			#set array element 1 for Zr
						cationPosition[m,3]=1
				
					if str(posPeriod[i,j,k,ii,0]) in cationList:
						#create neighbour list consisting of the 26 surround cells and the ijk cell
						neighbourList=np.concatenate((posPeriod[i,j,k,0:binCounter[i,j,k]],
													  posPeriod[i+1,j,k,0:binCounter[i+1,j,k]],posPeriod[i-1,j,k,0:binCounter[i-1,j,k]],
													  posPeriod[i,j+1,k,0:binCounter[i,j+1,k]],posPeriod[i,j-1,k,0:binCounter[i,j-1,k]],
													  posPeriod[i+1,j+1,k,0:binCounter[i+1,j+1,k]],posPeriod[i+1,j-1,k,0:binCounter[i+1,j-1,k]],
													  posPeriod[i-1,j+1,k,0:binCounter[i-1,j+1,k]],posPeriod[i-1,j-1,k,0:binCounter[i-1,j-1,k]],

													  posPeriod[i,j,k+1,0:binCounter[i,j,k+1]],
													  posPeriod[i+1,j,k+1,0:binCounter[i+1,j,k+1]],posPeriod[i-1,j,k+1,0:binCounter[i-1,j,k+1]],
													  posPeriod[i,j+1,k+1,0:binCounter[i,j+1,k+1]],posPeriod[i,j-1,k+1,0:binCounter[i,j-1,k+1]],
													  posPeriod[i+1,j+1,k+1,0:binCounter[i+1,j+1,k+1]],posPeriod[i+1,j-1,k+1,0:binCounter[i+1,j-1,k+1]],
													  posPeriod[i-1,j+1,k+1,0:binCounter[i-1,j+1,k+1]],posPeriod[i-1,j-1,k+1,0:binCounter[i-1,j-1,k+1]],

													  posPeriod[i,j,k-1,0:binCounter[i,j,k-1]],
													  posPeriod[i+1,j,k-1,0:binCounter[i+1,j,k-1]],posPeriod[i-1,j,k-1,0:binCounter[i-1,j,k-1]],
													  posPeriod[i,j+1,k-1,0:binCounter[i,j+1,k-1]],posPeriod[i,j-1,k-1,0:binCounter[i,j-1,k-1]],
													  posPeriod[i+1,j+1,k-1,0:binCounter[i+1,j+1,k-1]],posPeriod[i+1,j-1,k-1,0:binCounter[i+1,j-1,k-1]],
													  posPeriod[i-1,j+1,k-1,0:binCounter[i-1,j+1,k-1]],posPeriod[i-1,j-1,k-1,0:binCounter[i-1,j-1,k-1]],
													  ),axis=0)
						
						polarization[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4])			#define position of B cation
						tetFactor[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4])			
						unitVolume[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4])			
						strain[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4])	
						octVector[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4]); octAngle[m,0:3]=np.copy(posPeriod[i,j,k,ii,1:4])	

						PbCount=0; OCount=0; PbsCount=0; OsCount=0; TisCount=0; ZrsCount=0	#reset counters
						for jj in xrange(len(neighbourList)):

							if str(neighbourList[jj,0]) in ['Pb','O','Pbs','Os','Tis','Zrs','Mgs','Nbs']: 
								#calculate distance between B cation and all other ions
								distance=np.sqrt((float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))**2+(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))**2+(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))**2)    
								
								if distance<=5     and str(neighbourList[jj,0])=='Pb' :										#weight for Pb = 1/8
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))/8*PbCharge		#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))/8*PbCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))/8*PbCharge
									PbCount+=1 
									PbList[PbCount-1]=np.array([neighbourList[jj,1],neighbourList[jj,2],neighbourList[jj,3]])
								elif distance<=5   and str(neighbourList[jj,0])=='Pbs' :									#weight for Pb = 1/8
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))/8*PbsCharge		#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))/8*PbsCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))/8*PbsCharge
									PbsCount+=1	
								elif distance<=3   and str(neighbourList[jj,0])=='O' :										#weight for O = 1/2
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))/2*OCharge			#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))/2*OCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))/2*OCharge	
									OCount+=1	
									OList[OCount-1]=np.array([neighbourList[jj,1],neighbourList[jj,2],neighbourList[jj,3]])
								elif distance<=3   and str(neighbourList[jj,0])=='Os' :									#weight for O = 1/2
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))/2*OsCharge		#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))/2*OsCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))/2*OsCharge	
									OsCount+=1
								elif distance<=1.5 and str(neighbourList[jj,0])=='Tis' :									#weight for Ti = 1
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))*TisCharge			#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))*TisCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))*TisCharge	
									TisCount+=1
								elif distance<=1.5 and str(neighbourList[jj,0])=='Zrs' :									#weight for Zr = 1
									polarization[m,3]+=(float(neighbourList[jj,1])-float(posPeriod[i,j,k,ii,1]))*ZrsCharge			#add up polarization
									polarization[m,4]+=(float(neighbourList[jj,2])-float(posPeriod[i,j,k,ii,2]))*ZrsCharge
									polarization[m,5]+=(float(neighbourList[jj,3])-float(posPeriod[i,j,k,ii,3]))*ZrsCharge	
								 	ZrsCount+=1	 
								 							
						#calculate local lattice parameters and angles from PbList:
						aParameter=calcLattice(PbList,'x'); bParameter=calcLattice(PbList,'y'); cParameter=calcLattice(PbList,'z')
						unitVolume[m,3]=aParameter; unitVolume[m,4]=bParameter; unitVolume[m,5]=cParameter; unitVolume[m,6]=aParameter*bParameter*cParameter
						
						#calculate tetragonality factor
						tetFactor[m,3]=np.max(np.array([aParameter,bParameter,cParameter]))/np.min(np.array([aParameter,bParameter,cParameter]))
						
						#divide by local cell volume
						polarization[m,3]/=(aParameter*bParameter*cParameter); polarization[m,4]/=(aParameter*bParameter*cParameter); polarization[m,5]/=(aParameter*bParameter*cParameter)

						#calculate Octahedral Vector and Octahedral Angle		
						octVector[m,3:6]=calcOct(OList)[0]; octAngle[m,3:6]=calcOct(OList)[1]; 

						#calculate strain of each local cell
						strain[m,3]=(aParameter-(averageCellVectorA/superCell[0]))/(averageCellVectorA/superCell[0])*100
						strain[m,4]=(bParameter-(averageCellVectorB/superCell[1]))/(averageCellVectorB/superCell[1])*100
						strain[m,5]=(cParameter-(averageCellVectorC/superCell[2]))/(averageCellVectorC/superCell[2])*100
						strain[m,6]=((aParameter-(averageCellVectorA/superCell[0]))/(averageCellVectorB/superCell[1])+(bParameter-(averageCellVectorB/superCell[1]))/(averageCellVectorA/superCell[0]))*100*0.5
						strain[m,7]=((aParameter-(averageCellVectorA/superCell[0]))/(averageCellVectorC/superCell[2])+(cParameter-(averageCellVectorC/superCell[2]))/(averageCellVectorA/superCell[0]))*100*0.5
						strain[m,8]=((bParameter-(averageCellVectorB/superCell[1]))/(averageCellVectorC/superCell[2])+(cParameter-(averageCellVectorC/superCell[2]))/(averageCellVectorB/superCell[1]))*100*0.5
					
						#calculate strain of the whole cell TODO: Move out of loop
						cellStrain[0,0]=((cellVectorA[0]-averageCellVectorA)/(averageCellVectorA))*100
						cellStrain[0,1]=((cellVectorB[1]-averageCellVectorA)/(averageCellVectorB))*100
						cellStrain[0,2]=((cellVectorC[2]-averageCellVectorC)/(averageCellVectorC))*100
						cellStrain[0,3]=((cellVectorA[0]-averageCellVectorA)/averageCellVectorB+(cellVectorB[1]-averageCellVectorB)/averageCellVectorA)*100*0.5
						cellStrain[0,4]=((cellVectorA[0]-averageCellVectorA)/averageCellVectorC+(cellVectorC[2]-averageCellVectorC)/averageCellVectorA)*100*0.5
						cellStrain[0,5]=((cellVectorB[1]-averageCellVectorB)/averageCellVectorC+(cellVectorC[2]-averageCellVectorC)/averageCellVectorB)*100*0.5

						m+=1									#stop loop after all B cations have been found
						if m==nCell:
							break
						
						#check if 8 Pb and 6 O to make sure polarization is defined:
						if (PbCount+PbsCount+ OCount+ OsCount+TisCount+ZrsCount)!=29:
							polarization[m,3]=0; polarization[m,4]=0; polarization[m,5]=0;
							tetFactor[m,3]=1
							strain[m,3]=0; strain[m,4]=0; strain[m,5]=0;
							print "ERROR => No definied polarization at cell: " , polarization[m-1,0:3]
							
						#update progress bar
						if str(multiprocessing.current_process())=='<Process(PoolWorker-1, started daemon)>':	
							if m%(nCell*(((movie[1]+movie[2]-movie[0])/movie[2])/args.coreNumber)/toolbar_widh)==0:
								sys.stdout.write("-"); sys.stdout.flush()
	
	polarization[:,3:6]*=float(1.602177*10**(-19)*float(10**6))/(10**(-16))				#rescale polarization to microC/cm^2
	
	return polarization,tetFactor,unitVolume, strain,cellStrain,octVector,octAngle,cationPosition

def binPol(polarization,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor):
	"Divide system into cells and average polarization over those cells"
	#get size of input array (4 or 6)
	polArraySize=polarization.shape[1]
	#calculate bin sizes
	nCell=superCell[0]*superCell[1]*superCell[2]						#number of cells

	sizeA=np.linalg.norm(cellVectorA)/superCell[0]*binFactor[0]; sizeB=np.linalg.norm(cellVectorB)/superCell[1]*binFactor[1]; sizeC=np.linalg.norm(cellVectorC)/superCell[2]*binFactor[2]
	polPeriod=np.zeros(shape=(superCell[0]/binFactor[0],superCell[1]/binFactor[1],superCell[2]/binFactor[2],polArraySize),dtype=float)

	for m in xrange(nCell):												#loop over all polarization values

		polPeriod[math.floor(polarization[m,0]/sizeA),math.floor(polarization[m,1]/sizeB),math.floor(polarization[m,2]/sizeC),0:3]=[math.floor(polarization[m,0]/sizeA)*sizeA+2,math.floor(polarization[m,1]/sizeB)*sizeB+2, math.floor(polarization[m,2]/sizeC)*sizeC+2]

		polPeriod[math.floor(polarization[m,0]/sizeA),math.floor(polarization[m,1]/sizeB),math.floor(polarization[m,2]/sizeC),3:polArraySize]=polarization[m,3:polArraySize]

	polPeriod=polPeriod/(binFactor[0]*binFactor[1]*binFactor[2])

	if averageDirection!='N':
		directions={'x':0, 'X':0, 'y':1, 'Y':1, 'z':2, 'Z':2}							#dictionary for directions
		averageDirection=directions[averageDirection]
		#average along averageDirection
		polPeriod=np.mean(polPeriod,axis=averageDirection)#,weights=polPeriod.astype(bool))
		polAverage=np.reshape(polPeriod,(nCell/(superCell[averageDirection]*binFactor[1]*binFactor[2]),polArraySize),order='F')
		polAverage=polAverage[np.nonzero(np.nanmean(polAverage,axis=1))]		#remove elements that are completly zero to get only the elements corresponding to the chosen cation
		
	else:
		polAverage=np.reshape(polPeriod,(nCell,polArraySize),order='F')
		polAverage=polAverage[np.nonzero(np.average(polAverage,axis=1))]		#remove elements that are completly zero to get only the elements corresponding to the chosen cation

	if secondAverageDirection!=None:
		secondAverageDirection=directions[secondAverageDirection]
		#average along averageDirection
		polAverage=np.mean(polPeriod,axis=secondAverageDirection-1)

	return polAverage
	
def calcLattice(PbList,latticeDirection):
	"calculate local lattice parameters from 8 corner Pb atoms"
	directions={'x':0, 'X':0, 'y':1, 'Y':1, 'z':2, 'Z':2}									#dictionary for directions
	
	PbList=PbList[PbList[:,directions[latticeDirection]].argsort()]							#sort along direction
	topPb=np.mean(PbList[0:4,directions[latticeDirection]]); bottomPb=np.mean(PbList[4:8,directions[latticeDirection]])	#take 4 atoms each as top and bottom layer
	cellParameter=np.abs(topPb-bottomPb)										#calculate cell parameter difference between top and bottom layer
	
	return cellParameter
	
def calcOct(OList):
	"calculate the octahedron direction vector and angle"
	octDistance=np.empty([36,3])
	
	#calculate distances of Oxygen atoms
	k=0
	for i in xrange(6):
		for j in xrange(6):
			octDistance[k,0]=np.sqrt((OList[i,0]-OList[j,0])**2+(OList[i,1]-OList[j,1])**2+(OList[i,2]-OList[j,2])**2)
			octDistance[k,1]=i; octDistance[k,2]=j
			k+=1
			
	#select pair with largest distance as a first direction of the oxygen oxtahedron
	i=octDistance[np.argmax(octDistance[:,0])][1]
	j=octDistance[np.argmax(octDistance[:,0])][2]
	octVector=np.array([OList[i,0]-OList[j,0] , OList[i,1]-OList[j,1] ,OList[i,2]-OList[j,2]])
	
	octVector[0]=octVector[0]/np.linalg.norm(octVector); octVector[1]=octVector[1]/np.linalg.norm(octVector); octVector[2]=octVector[2]/np.linalg.norm(octVector)
	
	#fit a plane to the 4 planar Oxygen and to 2 out of plane Oxygen, then take the plane normal as Octahedron vector (biased towards out of plane octahedron (old direction)
	X=np.delete(OList,[i,j],axis=0)
	A=np.concatenate((2*X[:,0:3], [[1],[1],[1],[1]]),1)
	a,b,c,d = np.linalg.lstsq(A, np.append(octVector,(OList[i,2]+OList[j,2])/2))[0]																		#least square fit of octahedron

	newOctVector=np.copy(octVector)
	newOctVector[0]=a/np.linalg.norm(np.array([a,b,c])); newOctVector[1]=b/np.linalg.norm(np.array([a,b,c])); newOctVector[2]=c/np.linalg.norm(np.array([a,b,c]))	#normalize vector
	
	octVector=np.copy(newOctVector)*np.sign(np.max(newOctVector))																				#make vector facing in "positive" direction

	#calculate angles
	octAngle=np.empty(3)
	#alpha: angle with X axis, beta= angle with Y axis, gamma = angle with Z axis
	for i in range(3):
		if math.degrees(np.arccos(octVector[i]))>=60 and math.degrees(np.arccos(octVector[i]))<=150:
			octAngle[i]=math.degrees(np.arccos(octVector[i]))-90
		elif math.degrees(np.arccos(octVector[i]))>=150 :
			octAngle[i]=math.degrees(np.arccos(octVector[i]))-180
		else:
			octAngle[i]=math.degrees(np.arccos(octVector[i]))

	return octVector, octAngle

def printResult(result,name,frame):
	"print result array to file"
		
	np.set_printoptions(threshold=np.nan)											#print whole array

	#if only one cation selected change list accordingly:
	if cationSpecify=='Ti' or cationSpecify=='ti':
		cationString='Ti'
	elif cationSpecify=='Zr' or cationSpecify=='zr':
		cationString='Zr'
	elif cationSpecify=='default':
		cationString=''
	
	#print to file
	f = open(name+str(cationString)+'_'+str(frame)+'.dat', 'w')
	#print >>f, "#Timestep", frame
	print >>f, str(result).replace('[','').replace(']','')
	f.close()

def runMultiple(i):
	"run script for multiple timesteps"
	#initialize
	piezoList=np.empty([10]); piezoList[9]=i*0.0002	#timestep of 0.2 fs					#initialize list for polarization and strain with last entry as timestep

	m=range(startTimestep,totalTimestep*lengthTimestep,lengthTimestep).index(i)		#counter	for frames	
	historyFile='x'+str(range(totalTimestep).index(m)).zfill(5)
	
	cellVectorA, cellVectorB, cellVectorC=getCell(historyFile,i)							#get cell vectors

	position=getPos(historyFile,cellVectorA, cellVectorB, cellVectorC,totalAtom,i,superCell)        														#get positons
	polarization,tetFactor,unitVolume, strainFull,cellStrain,octVector,octAngle,cationPosition=calcPol(position,cellVectorA, cellVectorB, cellVectorC,i,superCell,cationSpecify,nCell) 		#calculate polarization
	polarization=binPol(polarization,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)						#bin polarization
	tetFactor=binPol(tetFactor,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)							#bin tetfactor
	unitVolume=binPol(unitVolume,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)						#bin unitVolume
	strain=binPol(strainFull,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)								#bin strain
	octVector=binPol(octVector, cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)							#bin octVector
	octAngle=binPol(octAngle, cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)							#bin octAngle
	cationPosition=binPol(cationPosition, cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)					#bin cationPosition
	
	printResult(polarization,'result', i); printResult(tetFactor,'tet', i);printResult(unitVolume,'vol', i); printResult(strain,'strain', i); printResult(octVector,'octVector', i); printResult(octAngle,'octAngle', i); 				#print to files 

	#create List of polarization and strain values for easy calculation of piezoelectric and dielectric tensors
	if os.path.isfile('OUTPUT') and outputTrue==1:					#if OUTPUT does not exist (or run didnt finish) use average of local strain instead of cellStrain
		piezoList[0:3]=np.nanmean(polarization,axis=0)[3:6]; piezoList[3:9]=cellStrain[0:6]
	else:
		piezoList[0:3]=np.nanmean(polarization,axis=0)[3:6]; piezoList[3:9]=np.nanmean(strainFull,axis=0)[3:9]

	return piezoList, polarization, tetFactor, unitVolume, strain, octVector, octAngle, cationPosition
		
def runSingle():
	#single timestep
	"get Cell Vectors  HISTORY file"
	historyFile=args.fileName
	inputFile=open(historyFile,'r')
        line=inputFile.readline()  
        
        while(line!=""): 
        	if string.split(line)[0][:8]=='timestep' and int(string.split(line)[1])==frame:
			line=inputFile.readline()
			cellVectorA=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])
			line=inputFile.readline()  
			cellVectorB=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])
			line=inputFile.readline()  
			cellVectorC=np.array([float(string.split(line)[0]),float(string.split(line)[1]),float(string.split(line)[2])])		
			break
		line=inputFile.readline()  						

	inputFile.close()

	atomList=['Pb', 'Ti', 'Zr','Mg','Nb', 'O','Pbs', 'Tis', 'Zrs','Os','Mgs','Nbs']
        position=np.empty(shape=(totalAtom,4),dtype=object)
        i=0                                                                    		 				#counter for position
        
        inputFile=open(historyFile,'r')
        line=inputFile.readline()  
	
        while(line!=""): 
        	if string.split(line)[0][:8]=='timestep' and int(string.split(line)[1])==frame:
        		line=inputFile.readline()
        		while(line!="") and i<totalAtom:
		                if string.split(line)[0][:3] in atomList:                       
		                        position[i,0]=str(string.split(line)[0][:3])             			    #get atom name
		                        
		                        line=inputFile.readline()                                      			    #advance a line
		                        
		                        position[i,1]=float(string.split(line)[0])%cellVectorA[0]                       #get x pos
		                        position[i,2]=float(string.split(line)[1])%cellVectorB[1]                       #get x pos
		                        position[i,3]=float(string.split(line)[2])%cellVectorC[2]                       #get x pos
		                        i+=1

		                line=inputFile.readline()
	        line=inputFile.readline() 	
        inputFile.close()

	print "Cell Dimensions: " ,'\n', cellVectorA,'\n', cellVectorB,'\n',cellVectorC,'\n'
	
	#get data, bin data, average data
	polarization,tetFactor,unitVolume,strain,cellStrain,octVector,octAngle,cationPosition=calcPol(position,cellVectorA, cellVectorB, cellVectorC,frame,superCell,cationSpecify,nCell) 
	polarization=binPol(polarization,cellVectorA, cellVectorB, cellVectorC, superCell,averageDirection,secondAverageDirection,binFactor)
	tetFactor=binPol(tetFactor,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)
	unitVolume=binPol(unitVolume,cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)
	octVector=binPol(octVector, cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)
	octAngle=binPol(octAngle, cellVectorA, cellVectorB, cellVectorC,superCell,averageDirection,secondAverageDirection,binFactor)
	bulkPol=np.nanmean(polarization,axis=0);bulkTet=np.nanmean(tetFactor,axis=0); bulkStrain=np.nanmean(strain,axis=0); 	
	bulkOctVector=(np.nanmean(octVector,axis=0)); bulkOctAngle=np.nanmean(octAngle,axis=0); 
	
	print "Bulk Polarization: ", bulkPol[3:6]
	print "Bulk Tetragonality: ", bulkTet[3:4]
	print "Bulk Strain: ", bulkStrain[3:9]
	print "Bulk Octahedron Vector (squared): ", bulkOctVector[3:6], "with angles: ", bulkOctAngle[3:6]
	
	printResult(polarization,'result', frame); printResult(tetFactor,'tet', frame); printResult(strain,'strain',frame); printResult(octVector,'octVector', i); printResult(octAngle,'octAngle', i); 	
	
###########################################################################################################################
#RUN MAIN ROUTINE

#get input from command line arguments (argparse on top of script)
frame=args.frame
averageDirection=args.averageDirection
secondAverageDirection=args.secondAverageDirection
timeAverage=args.timeAverage
inputFile=args.fileName
movie=args.movie
densvar=args.densvar
binFactor=args.binFactor
cationSpecify=args.cationSpecify
superCell=args.dim
correlationCalc=args.correlation
distributionCalc=args.distribution

nCell=superCell[0]*superCell[1]*superCell[2]		#number of cells
totalAtom, totalTimestep=getInfo(inputFile)			#number of atoms (from File)

#print info
print 'Total number of atoms: ', totalAtom/2
print 'Total number of cells: ', nCell
print "Supercell dimensions: ", superCell[0], superCell[1], superCell[2]
print 'Number of timesteps: ', totalTimestep

#loop over timesteps
if movie!=None and frame!=1:	
	print "You can't select TIMESTEP option together with MOVIE option. Choose either one!"
	sys.exit()
#loop over multiple timesteps
elif movie!=None:	
	
	#initialize
	piezoTensor=np.empty([3,6]); dielectricTensor=np.empty([1,3])	
	cellVectorA, cellVectorB, cellVectorC=getCell(inputFile,movie[1])					
	averageCellVectorA, averageCellVectorB, averageCellVectorC,outputTrue=getAverageCellVectors('OUTPUT',cellVectorA, cellVectorB, cellVectorC)	
	
	#split HISTORY into smaller parts (one per timestep) and get length of timestep
	print "\nSplit HISTORY file"
	os.system('tail -n +3 '+str(inputFile)+' | split -a 5 -d -l '+str(totalAtom*2+4)) 			
	f = open('x00000', 'r'); line=f.readline(); startTimestep=int(string.split(line)[1]); f.close()				#get starting timestep
	f = open('x00001', 'r'); line=f.readline(); lengthTimestep=int(string.split(line)[1])-startTimestep; f.close()  #get length of timesteps
		
	#calculate piezoList in paralell
	
	# setup progresssbar
	print "Calculate polarization:"
	toolbar_widh = 40
	sys.stdout.write("[%s]" % (" " * toolbar_widh)); sys.stdout.flush(); sys.stdout.write("\b" * (toolbar_widh+1)) 
		
	pool=Pool(args.coreNumber)
	resultArray=np.array(pool.map(runMultiple, range(movie[0],movie[1]+movie[2],movie[2])))			#get results
	pool.close(); pool.join(); 																	#join results 
	
	#split result array and sort arrays after timestep
	piezoList=np.vstack(resultArray[:,0])
	timeIndex=piezoList[:,9].argsort()	
	
	piezoList=np.vstack(resultArray[:,0])[timeIndex]
	polarization=np.vstack(resultArray[:,1]);	

	polarization=np.reshape(polarization,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(polarization)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),6))[timeIndex]
	tetFactor=np.vstack(resultArray[:,2]);		tetFactor=np.reshape(tetFactor,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(tetFactor)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),4))[timeIndex]
	unitVolume=np.vstack(resultArray[:,3]);	unitVolume=np.reshape(unitVolume,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(unitVolume)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),7))[timeIndex]
	strain=np.vstack(resultArray[:,4]);			strain=np.reshape(strain,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(strain)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),9),order='F')[timeIndex]
	octVector=np.vstack(resultArray[:,5]);		octVector=np.reshape(octVector,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(octVector)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),6))[timeIndex]
	octAngle=np.vstack(resultArray[:,6]);		octAngle=np.reshape(octAngle,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(octAngle)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),6))[timeIndex]
	cationPosition=np.vstack(resultArray[:,7]);	cationPosition=np.reshape(cationPosition,((movie[1]+movie[2]-movie[0])/movie[2],np.shape(cationPosition)[0]/((movie[1]+movie[2]-movie[0])/movie[2]),4))[timeIndex]
	cationPosition=np.nanmean(cationPosition,axis=0)		#time average cationPostion (Cation index (0 for Ti and 1 for Zr) should not change with time. This step is done to simplify the array (Ncell,4)
	
	#calculate covariance (piezoelectric coefficient)
	piezoTensor[0,0]=np.cov(piezoList[:,0],piezoList[:,3],rowvar=0)[0,1]; piezoTensor[0,1]=np.cov(piezoList[:,0],piezoList[:,4],rowvar=0)[0,1]; piezoTensor[0,2]=np.cov(piezoList[:,0],piezoList[:,5],rowvar=0)[0,1]
	piezoTensor[1,0]=np.cov(piezoList[:,1],piezoList[:,3],rowvar=0)[0,1]; piezoTensor[1,1]=np.cov(piezoList[:,1],piezoList[:,4],rowvar=0)[0,1]; piezoTensor[1,2]=np.cov(piezoList[:,1],piezoList[:,5],rowvar=0)[0,1]
	piezoTensor[2,0]=np.cov(piezoList[:,2],piezoList[:,3],rowvar=0)[0,1]; piezoTensor[2,1]=np.cov(piezoList[:,2],piezoList[:,4],rowvar=0)[0,1]; piezoTensor[2,2]=np.cov(piezoList[:,2],piezoList[:,5],rowvar=0)[0,1]
	piezoTensor[0,3]=np.cov(piezoList[:,0],piezoList[:,6],rowvar=0)[0,1]; piezoTensor[0,4]=np.cov(piezoList[:,0],piezoList[:,7],rowvar=0)[0,1]; piezoTensor[0,5]=np.cov(piezoList[:,0],piezoList[:,8],rowvar=0)[0,1]
	piezoTensor[1,3]=np.cov(piezoList[:,1],piezoList[:,6],rowvar=0)[0,1]; piezoTensor[1,4]=np.cov(piezoList[:,1],piezoList[:,7],rowvar=0)[0,1]; piezoTensor[1,5]=np.cov(piezoList[:,1],piezoList[:,8],rowvar=0)[0,1]
	piezoTensor[2,3]=np.cov(piezoList[:,2],piezoList[:,6],rowvar=0)[0,1]; piezoTensor[2,4]=np.cov(piezoList[:,2],piezoList[:,7],rowvar=0)[0,1]; piezoTensor[2,5]=np.cov(piezoList[:,2],piezoList[:,8],rowvar=0)[0,1]

	#get temperature to calculate Piezo and Dielectric Tensor
	foundLine=0
	if os.path.isfile('OUTPUT'):
		searchFile=open('OUTPUT','r')
		for line in searchFile:
			if line.startswith(' run terminated after'):
				#skip lines
				line=searchFile.next();line=searchFile.next();line=searchFile.next();line=searchFile.next();line=searchFile.next();line=searchFile.next();
				line=searchFile.next();line=searchFile.next();line=searchFile.next();line=searchFile.next();
				simulationTemp=float(line.split()[2])
				foundLine=1
				break
		
	elif os.path.isfile('OUTPUT') and foundLine==0:
		print "Can't access temperature from OUTPUT file!"
		print "--> temperature = 200 K is used to calculate Piezo and Dielectric Tensors" ,'\n'		
		simulationTemp=200 								#set temperature to200 K
	else:
		simulationTemp=200
		
	piezoTensor*=float(((averageCellVectorA*averageCellVectorB*averageCellVectorC*10**-30)*10**-4)/(simulationTemp*1.380*10**-23)*10**12)			#normalize to pC/N
						
	#calculate dielectric Tensor
	dielectricTensor[0,0]=np.mean(piezoList[:,0]**2,axis=0)-np.mean(piezoList[:,0],axis=0)**2
	dielectricTensor[0,1]=np.mean(piezoList[:,1]**2,axis=0)-np.mean(piezoList[:,1],axis=0)**2
	dielectricTensor[0,2]=np.mean(piezoList[:,2]**2,axis=0)-np.mean(piezoList[:,2],axis=0)**2
	dielectricTensor*=float(((averageCellVectorA*averageCellVectorB*averageCellVectorC*10**-30)*10**-4)/(8.854*10**-12*simulationTemp*1.380*10**-23))	#normalize 
	
	#PRINT
	#if only one cation selected change list accordingly:
	if cationSpecify=='Ti' or cationSpecify=='ti':
		cationString='Ti'
	elif cationSpecify=='Zr' or cationSpecify=='zr':
		cationString='Zr'
	elif cationSpecify=='default':
		cationString=''
	
	np.set_printoptions(threshold=np.nan)						#print whole array
	
	#print convergence of d33 and e33 to conv_piezo.dat to check convergence of these values
	f = open('conv_piezo'+str(cationString)+'.dat', 'w')
	print>>f, '#timestep , d33 , e33'
		
	for i in xrange(np.shape(piezoList)[0]):
		print >>f, str(i), str(np.cov(piezoList[0:i+1,2],piezoList[0:i+1,5],rowvar=0)[0,1]*float(((averageCellVectorA*averageCellVectorB*averageCellVectorC*10**-30)*10**-4)/(simulationTemp*1.380*10**-23)*10**12)), str((np.mean(piezoList[0:i+1,2]**2,axis=0)-np.mean(piezoList[0:i+1,2],axis=0)**2)*float(((averageCellVectorA*averageCellVectorB*averageCellVectorC*10**-30)*10**-4)/(8.854*10**-12*simulationTemp*1.380*10**-23)))
	f.close()	

	#print polarization and strain at each timestep to file. 
	f = open('polarization'+str(cationString)+'.dat', 'w')
	print>>f, '#timestep, Px, Py, Pz, Exx, Eyy, Ezz, Exy, Exz, Eyz'
	print>>f, str(piezoList[:,0:10]).replace('[','').replace(']','')
	f.close()

	#print piezoTensor to file
	f = open('piezo'+str(cationString)+'.dat', 'w')
	print>>f, '#dxx, dyy, dzz, dxy, dxz, dyz'
	print>>f, str(piezoTensor).replace('[','').replace(']','')
	f.close()

	#print dielectricTensor to file
	f = open('dielectric'+str(cationString)+'.dat', 'w')
	print>>f, '#Dxx, Dyy, Dzz'
	print>>f, str(dielectricTensor).replace('[','').replace(']','')
	f.close()

	#average over timesteps
	if timeAverage==True:
		print "\n\nCalculating time average over",((movie[2]+movie[1]-movie[0])/movie[2]), "timesteps\n"	

		bulkPol=np.nanmean(piezoList,axis=0); bulkPol=np.append(bulkPol, np.sqrt(bulkPol[0]**2+bulkPol[1]**2+bulkPol[2]**2))	#calculate bulk polarization as average
		bulkStrain=np.nanmean(piezoList,axis=0)																	#calculate bulk strain as average
				
		print "Bulk Polarization:      ", bulkPol[0:3], np.sqrt(bulkPol[0]**2+bulkPol[1]**2+bulkPol[2]**2); printResult(bulkPol[0:3],'result','bulk')
		print "Bulk Strain:            ", bulkStrain[3:9] ,'\n' ; printResult(bulkPol[3:9],'strain', 'bulk')

		print "Bulk octahedral angles: ", np.nanmean(np.nanmean(octAngle,axis=0),axis=0)[3:6],'and rms: ', np.nanmean(np.nanmean(np.sqrt(octAngle**2),axis=0),axis=0)[3:6] 
		print "Bulk octahedral vector: ", np.nanmean(np.nanmean(octVector,axis=0),axis=0)[3:6] ,'\n'

		print "\nPiezoelectric Tensor (calculated at ", simulationTemp, "K): \n ", str(piezoTensor).replace('[','').replace(']','') ,'\n'	
		print "Dielectric Tensor (calculated at ", simulationTemp, "K): \n ", str(dielectricTensor).replace('[','').replace(']','') ,'\n'

		
		#PRINT
		printResult(np.nanmean(polarization,axis=0),'result','averaged')
		printResult(np.nanmean(strain,axis=0),'strain','averaged')
		printResult(np.nanmean(tetFactor,axis=0),'tet','averaged'); printResult(np.nanmean(np.nanmean(tetFactor,axis=0),axis=0),'tet','bulk')
		printResult(np.nanmean(unitVolume,axis=0),'vol','averaged'); printResult(np.nanmean(np.nanmean(unitVolume,axis=0),axis=0),'vol','bulk')
		printResult(np.nanmean(octVector,axis=0),'octVector','averaged'); printResult(np.nanmean(np.nanmean(octVector,axis=0),axis=0)[3:6],'octVector','bulk')
		printResult(np.nanmean(octAngle,axis=0),'octAngle','averaged'); printResult(np.nanmean(np.nanmean(octAngle,axis=0),axis=0)[3:6],'octAngle','bulk')
		printResult(np.nanmean(np.nanmean(np.sqrt(octAngle**2),axis=0),axis=0)[3:6],'octAngle','bulkrms')
	
	#calculate polarization distribution
	if distributionCalc!=None and timeAverage==True:
		print "Calculate polarization distributions...\n"
		distPolX=np.histogram(np.nanmean(polarization,axis=0)[:,3],bins=100)[1][np.newaxis][:,0:100]
		distPolX=np.append(distPolX,np.array(np.histogram(np.nanmean(polarization,axis=0)[:,3],bins=100)[0])[np.newaxis],axis=0)
		distPolY=np.histogram(np.nanmean(polarization,axis=0)[:,4],bins=100)[1][np.newaxis][:,0:100]
		distPolY=np.append(distPolY,np.array(np.histogram(np.nanmean(polarization,axis=0)[:,4],bins=100)[0])[np.newaxis],axis=0)	
		distPolZ=np.histogram(np.nanmean(polarization,axis=0)[:,5],bins=100)[1][np.newaxis][:,0:100]
		distPolZ=np.append(distPolZ,np.array(np.histogram(np.nanmean(polarization,axis=0)[:,5],bins=100)[0])[np.newaxis],axis=0)	
		distPolTot=np.histogram(np.sqrt(np.nanmean(polarization,axis=0)[:,3]**2+np.nanmean(polarization,axis=0)[:,4]**2+np.nanmean(polarization,axis=0)[:,5]**2),bins=100)[1][np.newaxis][:,0:100]
		distPolTot=np.append(distPolTot,np.array(np.histogram(np.sqrt(np.nanmean(polarization,axis=0)[:,3]**2+np.nanmean(polarization,axis=0)[:,4]**2+np.nanmean(polarization,axis=0)[:,5]**2),bins=100)[0])[np.newaxis],axis=0)	

		printResult(distPolX.T,'distX','averaged'); printResult(distPolY.T,'distY','averaged'); printResult(distPolZ.T,'distZ','averaged'); printResult(distPolTot.T,'distTot','averaged');

		#print out xyz files sorted after polarization
		id=np.transpose(np.where((np.nanmean(polarization,axis=0)[:,5]> distributionCalc[0] )*(np.nanmean(polarization,axis=0)[:,5]< distributionCalc[1] )))
		distPolPos=np.nanmean(polarization,axis=0)[id,0:3][:,0]
	
		x = open('Pos_'+str(distributionCalc[0])+'_'+str(distributionCalc[1])+'.xyz', 'w')
		print >>x,str(np.size(distPolPos)/3)
		for n in xrange(np.size(distPolPos)/3):
			print >>x , 'Fe', str(distPolPos[n]).replace('[','').replace(']','')
		x.close()	

		#calulcate nearest neighbours if --distribution switch is (float, float , 1):
		if distributionCalc[2]==1:
		
			#initialize periodic arrays
			cationPeriod=np.empty(shape=(superCell[0],superCell[1],superCell[2],4),dtype=object)
			cationPositionLargePeriod=np.empty(shape=(superCell[0]+2,superCell[1]+2,superCell[2]+2,4),dtype=object)
			cationPositionLarge=np.empty(((superCell[0]+2)*(superCell[1]+2)*(superCell[2]+2),4),dtype=object)
			binCounter=np.ones(shape=(superCell[0],superCell[1],superCell[2]),dtype=int)					#Counter for number of atoms in each bin

			#make cationPosition 3d periodic array
			m=0										
			for s in xrange(superCell[2]):
				for r in xrange(superCell[1]):
					for q in xrange(superCell[0]):
						cationPeriod[q,r,s,0]=cationPosition[m,3]
						cationPeriod[q,r,s,1:4]=cationPosition[m,0:3]
						m+=1

			#extend cationPosition array by 2 unit cells in each direction (to avoid PBC)			
			cationPositionLargePeriod[1:superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1]=np.copy(cationPeriod) 					#copy position array
			#along X
			cationPositionLargePeriod[0,1:superCell[1]+1,1:superCell[2]+1,1]=cationPositionLargePeriod[superCell[0],1:superCell[1]+1,1:superCell[2]+1,1]-cellVectorA[0] 
			cationPositionLargePeriod[0,1:superCell[1]+1,1:superCell[2]+1,2]=cationPositionLargePeriod[superCell[0],1:superCell[1]+1,1:superCell[2]+1,2]
			cationPositionLargePeriod[0,1:superCell[1]+1,1:superCell[2]+1,3]=cationPositionLargePeriod[superCell[0],1:superCell[1]+1,1:superCell[2]+1,3]
			cationPositionLargePeriod[0,1:superCell[1]+1,1:superCell[2]+1,0]=cationPositionLargePeriod[superCell[0],1:superCell[1]+1,1:superCell[2]+1,0]  
			cationPositionLargePeriod[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1,1]=cationPositionLargePeriod[1,1:superCell[1]+1,1:superCell[2]+1,1]+cellVectorA[0] 
			cationPositionLargePeriod[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1,2]=cationPositionLargePeriod[1,1:superCell[1]+1,1:superCell[2]+1,2]
			cationPositionLargePeriod[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1,3]=cationPositionLargePeriod[1,1:superCell[1]+1,1:superCell[2]+1,3] 
			cationPositionLargePeriod[superCell[0]+1,1:superCell[1]+1,1:superCell[2]+1,0]=cationPositionLargePeriod[1,1:superCell[1]+1,1:superCell[2]+1,0] 
			#along Y	
			cationPositionLargePeriod[:,0,1:superCell[2]+1,1]=cationPositionLargePeriod[:,superCell[1],1:superCell[2]+1,1]; 
			cationPositionLargePeriod[:,0,1:superCell[2]+1,2]=cationPositionLargePeriod[:,superCell[1],1:superCell[2]+1,2]-cellVectorB[1] ;  
			cationPositionLargePeriod[:,0,1:superCell[2]+1,3]=cationPositionLargePeriod[:,superCell[1],1:superCell[2]+1,3];  
			cationPositionLargePeriod[:,0,1:superCell[2]+1,0]=cationPositionLargePeriod[:,superCell[1],1:superCell[2]+1,0];    
			cationPositionLargePeriod[:,superCell[1]+1,1:superCell[2]+1,1]=cationPositionLargePeriod[:,1,1:superCell[2]+1,1]; 
			cationPositionLargePeriod[:,superCell[1]+1,1:superCell[2]+1,2]=cationPositionLargePeriod[:,1,1:superCell[2]+1,2]+cellVectorB[1] ; 
			cationPositionLargePeriod[:,superCell[1]+1,1:superCell[2]+1,3]=cationPositionLargePeriod[:,1,1:superCell[2]+1,3]; 
			cationPositionLargePeriod[:,superCell[1]+1,1:superCell[2]+1,0]=cationPositionLargePeriod[:,1,1:superCell[2]+1,0]; 
			#along Z
			cationPositionLargePeriod[:,:,0,1]=cationPositionLargePeriod[:,:,superCell[2],1];  
			cationPositionLargePeriod[:,:,0,2]=cationPositionLargePeriod[:,:,superCell[2],2];  
			cationPositionLargePeriod[:,:,0,3]=cationPositionLargePeriod[:,:,superCell[2],3]-cellVectorC[2] ;  
			cationPositionLargePeriod[:,:,0,0]=cationPositionLargePeriod[:,:,superCell[2],0];  
			cationPositionLargePeriod[:,:,superCell[2]+1,1]=cationPositionLargePeriod[:,:,1,1]; 
			cationPositionLargePeriod[:,:,superCell[2]+1,2]=cationPositionLargePeriod[:,:,1,2]; 
			cationPositionLargePeriod[:,:,superCell[2]+1,3]=cationPositionLargePeriod[:,:,1,3]+cellVectorC[2] ;   
			cationPositionLargePeriod[:,:,superCell[2]+1,0]=cationPositionLargePeriod[:,:,1,0]; 
			
			#flatten cationPositionLargePeriod 
			m=0
			for s in xrange(superCell[2]+2):
				for r in xrange(superCell[1]+2):
					for q in xrange(superCell[0]+2):
						cationPositionLarge[m,3]=cationPositionLargePeriod[q,r,s,0]
						cationPositionLarge[m,0:3]=cationPositionLargePeriod[q,r,s,1:4]
						m+=1
						
			#Build reference list with cation Name for polarization distributions
			cationReference=np.empty([np.size(distPolPos)/3,4],dtype='object')
			for i in xrange(np.size(distPolPos)/3):
				for j in xrange(nCell):
					if (cationPosition[j,0]==distPolPos[i,0]  and cationPosition[j,1]==distPolPos[i,1]  and cationPosition[j,2]==distPolPos[i,2]):
						cationReference[i]=cationPosition[j]

			#find Neighbours of same species in polarization distributions
			print "Finding nearest neighbours with Z polarization values from: " ,distributionCalc[0], 'to ', distributionCalc[1] 
			
			# setup progresss bar
			toolbar_widh = 40; sys.stdout.write("[%s]" % (" " * toolbar_widh)); sys.stdout.flush(); sys.stdout.write("\b" * (toolbar_widh+1)) 
		
			TiNeighbour=[]; ZrNeighbour=[]; TiNextNeighbour=[]; ZrNextNeighbour=[]
			for i in xrange(np.size(cationReference)/4):
				neighbourCounter=0; nextNeighbourCounter=0
				for j in xrange(np.size(cationPositionLarge)/4): 
					cationDistance=np.sqrt((cationReference[i,0]-cationPositionLarge[j,0])**2+(cationReference[i,1]-cationPositionLarge[j,1])**2+(cationReference[i,2]-cationPositionLarge[j,2])**2) #calculate distance
					
					#count nearest neighbours
					if cationDistance <= 4.9 and cationDistance > 0.1:  																						#avoid counting the center ion itself
						if cationReference[i,3]==cationPositionLarge[j,3]:
							neighbourCounter+=1
					#count next nearest neighbours	
					if cationDistance <= 6.2 and cationDistance > 0.1:  																							#avoid counting the center ion itself
						if cationReference[i,3]==cationPositionLarge[j,3]:
							nextNeighbourCounter+=1
						
				if cationReference[i,3]==0:
					TiNeighbour=np.append(TiNeighbour,neighbourCounter)
					TiNextNeighbour=np.append(TiNextNeighbour,nextNeighbourCounter)
				elif cationReference[i,3]==1:		
					ZrNeighbour=np.append(ZrNeighbour,neighbourCounter)
					ZrNextNeighbour=np.append(ZrNextNeighbour,nextNeighbourCounter)
					
				#update progressbar	
				if i%(np.shape(cationReference)[0]/toolbar_widh*2)==0:
					sys.stdout.write("-"); sys.stdout.flush()
		
			#print results
			print "Average number of nearest neighbours atoms withing polarization values: " ,distributionCalc[0], 'to ' ,distributionCalc[1] ,':'
			print "Ti: ", np.mean(TiNeighbour), "\tZr: ", np.mean(ZrNeighbour), "\tAll: ", np.mean(np.append(TiNeighbour,ZrNeighbour) )
			print "Ti: ", np.mean(TiNeighbour)/6*100, '(%)', "\tZr: ", np.mean(ZrNeighbour)/6*100,  '(%)',"\tAll: ", np.mean(np.append(TiNeighbour,ZrNeighbour) )/6*100 ,'(%)'
			print "Average number of next nearest neighbours atoms withing polarization values: " ,distributionCalc[0], 'to ' ,distributionCalc[1] ,':'
			print "Ti: ", np.mean(TiNextNeighbour), "\tZr: ", np.mean(ZrNextNeighbour), "\tAll: ", np.mean(np.append(TiNextNeighbour,ZrNextNeighbour) )
			print "Ti: ", np.mean(TiNextNeighbour)/18*100,  '(%)', "\tZr: ", np.mean(ZrNextNeighbour)/18*100,  '(%)',"\tAll: ", np.mean(np.append(TiNextNeighbour,ZrNextNeighbour) )/18*100, '(%)'
			#print to files
			f = open('neighbour_'+str(distributionCalc[0])+'_'+str(distributionCalc[1])+'.dat', 'w')
			print >>f, "Average number of nearest neighbours atoms withing polarization values: " ,distributionCalc[0], 'to ' ,distributionCalc[1] ,':'
			print >>f, "Ti: ", np.mean(TiNeighbour), "\tZr: ", np.mean(ZrNeighbour), "\tAll: ", np.mean(np.append(TiNeighbour,ZrNeighbour) )
			print >>f, "Ti: ", np.mean(TiNeighbour)/6*100, '(%)', "\tZr: ", np.mean(ZrNeighbour)/6*100,  '(%)',"\tAll: ", np.mean(np.append(TiNeighbour,ZrNeighbour) )/6*100, '(%)'
			print >>f, "Average number of next  nearest neighbours atoms withing polarization values: " ,distributionCalc[0], 'to ' ,distributionCalc[1] ,':'
			print >>f, "Ti: ", np.mean(TiNextNeighbour), "\tZr: ", np.mean(ZrNextNeighbour), "\tAll: ", np.mean(np.append(TiNextNeighbour,ZrNextNeighbour) )
			print >>f, "Ti: ", np.mean(TiNextNeighbour)/18*100, '(%)', "\tZr : ", np.mean(ZrNextNeighbour)/18*100, '(%)', "\tAll: ", np.mean(np.append(TiNextNeighbour,ZrNextNeighbour) )/18*100 ,'(%)'
			f.close()
			
			#for debugging
#			neighHist=np.histogram(TiNeighbour,bins=100)[1][np.newaxis][:,0:100]
#			neighHist=np.append(neighHist,np.array(np.histogram(TiNeighbour,bins=100)[0])[np.newaxis],axis=0)
#			ff = open('neighHist.dat', 'w')
#			for ll in xrange(np.size(neighHist)/2):
#				print >>ff, str(neighHist.T[ll]).replace('[','').replace(']','')
			

	#calculate spatial correlation
	if correlationCalc==True and timeAverage!=True:
		print "Select time averaging (--timeAverage) as option in order to calculate corrleations. Beware of averaging along space coordinates"
		
	elif	correlationCalc==True and timeAverage==True:
		print "Calculate spatial correlation:"

		averagedPos=np.nanmean(polarization,axis=0)[:,0:3]
		averagedPolarization=np.nanmean(polarization,axis=0)[:,0:6]; stdPolarization=np.std(averagedPolarization,axis=0)
		averagedOctVector=np.nanmean(octVector,axis=0)[:,0:6]; stdOctVector=np.std(averagedOctVector,axis=0)
		
		distanceArray=np.arange(0,cellVectorC[2]*2,0.1)
		correlationPolarization=np.zeros(np.size(distanceArray)); correlationOctVector=np.zeros(np.size(distanceArray));
		counterArray=np.zeros(np.size(distanceArray))

		# setup progresss bar
		toolbar_widh = 40
		sys.stdout.write("[%s]" % (" " * toolbar_widh)); sys.stdout.flush(); sys.stdout.write("\b" * (toolbar_widh+1)) 
		
		#slow double loop: TODO parallize this loop over i 
		for i in xrange(np.shape(averagedPos)[0]):
			for j in xrange(i,np.shape(averagedPos)[0]):

				k=np.round(np.sqrt((averagedPos[i,0]-averagedPos[j,0])**2+(averagedPos[i,1]-averagedPos[j,1])**2+(averagedPos[i,2]-averagedPos[j,2])**2),decimals=1)*10		#find index in distanceArray
				
				correlationPolarization[k]+=np.sqrt(((averagedPolarization[i,3]-stdPolarization[3])*(averagedPolarization[j,3]-stdPolarization[3]))**2+((averagedPolarization[i,4]-stdPolarization[4])*(averagedPolarization[j,4]-stdPolarization[4]))**2+((averagedPolarization[i,5]-stdPolarization[5])*(averagedPolarization[j,5]-stdPolarization[5]))**2)
				correlationOctVector[k]+=((averagedOctVector[i,3]-stdOctVector[3])*(averagedOctVector[j,3]-stdOctVector[3]))+((averagedOctVector[i,4]-stdOctVector[4])*(averagedOctVector[j,4]-stdOctVector[4]))+((averagedOctVector[i,5]-stdOctVector[5])*(averagedOctVector[j,5]-stdOctVector[5]))
				counterArray[k]+=1
				
			#update progressbar	
			if i%(np.shape(averagedPos)[0]/toolbar_widh)==0:
				sys.stdout.write("-"); sys.stdout.flush()
			
		#remove zero distance elements
		correlationPolarizationOld=np.copy(correlationPolarization); correlationOctVector=np.copy(correlationOctVector)
		correlationPolarization=correlationPolarization[np.nonzero(correlationPolarizationOld)]; 	correlationPolarization=np.delete(correlationPolarization,0)
		correlationOctVector=correlationOctVector[np.nonzero(correlationOctVector)]; 			correlationOctVector=np.delete(correlationOctVector,0)
		
		distanceArray=distanceArray[np.nonzero(correlationPolarizationOld)]; 				distanceArray=np.delete(distanceArray,0)
		counterArray=counterArray[np.nonzero(correlationPolarizationOld)]; 					counterArray=np.delete(counterArray,0)
		
		#normalizing...correlations
		correlationPolarization=correlationPolarization/counterArray; correlationOctVector=correlationOctVector/counterArray
		correlationPolarization=correlationPolarization/np.sum(correlationPolarization);	correlationOctVector=correlationOctVector/np.sum(correlationOctVector);	
	
		correlationPolarization=np.vstack((distanceArray,correlationPolarization)).T; correlationOctVector=np.vstack((distanceArray,correlationOctVector)).T; 

		#print to files
		np.set_printoptions(threshold=np.nan)		#print whole array
		f = open('correlation_polarization'+str(cationString)+'.dat', 'w')
		print>>f, "#distance, correlation"
		print>>f, str(correlationPolarization).replace('[','').replace(']','')
		f.close()
		
		f = open('correlation_octVector'+str(cationString)+'.dat', 'w')
		print>>f, "#distance, correlation"
		print>>f, str(correlationOctVector).replace('[','').replace(']','')
		f.close()

	#cleanup 
	os.system('rm x*')
	if os.path.isdir('data'): 
		os.system('mv *0.dat data')
	else:
		os.system('mkdir data'); os.system('mv *0.dat data')
	
#single timestep
else:	
	runSingle()
