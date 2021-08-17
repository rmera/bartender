#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

"""
Compare.py a scrit to analyze itp outputs from bartender

use:

/path/Compare.py gmx1.itp gmx2.itp ...  [-iplot]


It will produce boxplots comparing the parameters from each file starting from the second, to those from the first file.
with the optional flag -iplot, it will show the interactive plots. Otherwise the figures as just saved in the current dir.
"""


def Param():
    return {"data":[],"rmsd":-1,"atoms":[]}

class read:
    def __init__(self):
        self.r={"bonds":False,"angles":False,"reb":False,"improp":False,"dihe":False}
    def false(self):
        for i in self.r.keys():
            self.r[i]=False
    def w(self):
        for i in self.r.keys():
            if self.r[i]:
                return i
        



def fill_hooke(line,atoms):
        b=Param()
        l=line.split()
        for i in range(atoms):
            b["atoms"].append(int(l[i]))
        b["data"].append(float(l[-5]))
        b["data"].append(float(l[-4]))
        b["rmsd"]=float(l[-1])
        return b
            


"""
gmx_out_reader parses the data in a bartender itp file. The bartender part is important, because it expects the comments 
that Bartender adds. If you want to read a non-bartender itp (say, one of the existing, Groningen topologies for
comparison) you have to make the formats match. It takes into an argument a list, and adds one more element to that list.
The element added is the datosq dictionary defined in the first line of the function.

Each element of the dictionary is a list (one element per bond, or angle, etc. in the file) of dictionaries.

Each of these dictionaries represents a bond or an angle or dihedral, etc., and contains the following elements:

data:  a list of floats with the parameters for the element: eq distance or angles, force constant and (for simple dihedrals), periodicity, in that order
or the 6 C coefficients for a Ryckaert-Bellemans potential.

atoms: A list with the indexes (as read from the file, i.e. counting from one) of the atoms involved in the element

rmsd: a float, the rmsd for the Bartender fit for that element.

If you want to add more analyses, you need to build a function that can read gmx_out_reader output.
"""
def gmx_out_reader(filename,dql):
    datosq={"Reb":[],"Bonds":[],"Angles":[],"Simpled":[],"RB":[],"Improp":[]}
    dql.append(datosq)
    r=read()
    fin=open(filename,"r")
    first=False
    for i in fin:
        if "[bonds]" in i or "[constraints]" in i:
            r.false()
            r.r["bonds"]=True
            first=True
            continue
        if "[angles]" in i:
            r.false()
            r.r["angles"]=True
            first=True
            continue
        if "; Reb" in i:
            r.false()
            r.r["reb"]=True
            first=True
            continue
        if ";Improper" in i:
            r.false()
            r.r["improp"]=True
            first=True
            continue
        if "[dihedrals]" in i:
            r.false()
            r.r["dihe"]=True
            first=True
            continue
        if r.w()=="bonds":
            if i.startswith(";"):
                continue
            dql[-1]["Bonds"].append(fill_hooke(i,2))
        if r.w()=="bonds":
            if i.startswith(";"):
                continue
            dql[-1]["Bonds"].append(fill_hooke(i,2))
        if r.w()=="angles":
            if i.startswith(";"):
                continue
            dql[-1]["Angles"].append(fill_hooke(i,3))
        if r.w()=="reb":
            if i.startswith(";"):
                continue
            dql[-1]["Reb"].append(fill_hooke(i,3))
        if r.w()=="improp":
            if i.startswith(";"):
                continue
            dql[-1]["Improp"].append(fill_hooke(i,4))
        if r.w()=="dihe":
            if i.startswith(";;;"):
                l=i.lstrip(";;;").split()
                rb=Param()
                for j in range(4):
                    rb["atoms"].append(int(l[j]))
                rb["rmsd"]=float(l[-1])
                rb["data"].append(float(l[5]))
                rb["data"].append(float(l[6]))
                rb["data"].append(float(l[7]))
                rb["data"].append(float(l[8]))
                rb["data"].append(float(l[9]))
                rb["data"].append(float(l[10]))
                dql[-1]["RB"].append(rb)
                continue
            if i.startswith(";"):
                continue
            l=i.split()
            dih=Param()
            dih["atoms"]=[]
            for j in range(4):
                dih["atoms"].append(int(l[j]))
            dih["rmsd"]=float(l[-1])
            dih["data"]=[]
            dih["data"].append(float(l[-6]))
            dih["data"].append(float(l[-5]))
            dih["data"].append(int(l[-4])) # n
            dql[-1]["Simpled"].append(dih)

def proc_hooke(data):
    aggre=[[],[]]
    for j,w in enumerate(data):
        aggre[0].append(data[j]["data"][0])
    for j,w in enumerate(data):
        aggre[1].append(data[j]["data"][1])
    for i,v  in enumerate(aggre):
        aggre[i]=np.array(aggre[i])
    return aggre[0],aggre[1]


def proc_rb(data):
    aggre=[]
    for j,w in enumerate(data):
        for i,v in enumerate(w):
            aggre.append(data[j]["data"][i])
    aggre=np.array(aggre)
    return aggre





"""
CompareMethods takes a list produced by gmx_out_read and compares all elements (except for the first) to the first one.
It produces box plot for each type of bonded term (one plot for all the bonds, etc) where the deviations for each element of the term
(equilibrium distances or angles, and force constants) from the corresponding term in the first element are grouped over all bonds
and 3 box plots are produced for each term: One for the eq. distances of angles (one box per itp file analyzed), one for the force constants
and one for the RMSDs (the RMSDs plot contains one more box, because there the first itp, the reference, is also included).
Note that the periodicity of simple diedrals (Simpled) is not considered.
For R-B terms, 2 plots are produced: One for all the  C coefficients grouped together, and for the RMSD
"""
def CompareMethods(methods):
    #first we take all the data and subtract the equivalent for the reference, and square it
    ref=methods[0]
    for i,v in enumerate(methods[1:]):
        #this has got to be my spaghetti-est code ever ^_^
        for j in v.keys():
            for l,ww in enumerate(v[j]):
                if ref[j][l]["atoms"] != v[j][l]["atoms"]:
                    raise("parameters are not in the same order, quedo la mansaca!!") #this shouldn't happen. If it does, I have to rewrite this function.
                for k,val in enumerate(methods[i+1][j][l]["data"]):
                        d=methods[i+1][j][l]["data"][k]
                        refd=ref[j][l]["data"][k]
                        methods[i+1][j][l]["data"][k]=(d-refd)
    #each key will be a plot    
    fsize=(20,9)
    plots_hooke=["Reb","Bonds","Angles","Improp","Simpled"]
    for i in plots_hooke:
        eq=[]
        k=[]
        rmsd=[]
        if  len(methods[0][i])==0:
            continue
        for index,j in enumerate(methods):
            rmsdsmall=[]
            for l in j[i]:
                rmsdsmall.append(l["rmsd"])
            rmsd.append(np.array(rmsdsmall))
            if index==0:
                continue
            eqn,kn=proc_hooke(j[i])
            eq.append(eqn)
            k.append(kn)
            
        #create the plots.
        fig, ax = plt.subplots(3,1,figsize=fsize)
        fig.suptitle(i)
        ax[0].set_title("Equilibrium values")
        ax[0].boxplot(eq)

        ax[1].set_title("Force Constants")
        ax[1].boxplot(k)

        ax[2].set_title("RMSD")
        ax[2].boxplot(rmsd)
        plt.savefig(i+".png")
        if "-iplot" in sys.argv:
            plt.show()

### The RB plots

    CS=[]
    rmsd=[]
    i="RB"
    if  len(methods[0][i])==0:
        return
    for index,j in enumerate(methods):
        rmsdsmall=[]
        for l in j[i]:
            rmsdsmall.append(l["rmsd"])
        rmsd.append(np.array(rmsdsmall))
        if index==0:
            continue
        cs=proc_rb(j[i])
        CS.append(cs)  
    #create the plots.
    fig, ax = plt.subplots(2,1,figsize=fsize)
    fig.suptitle(i)
    ax[0].set_title("Coefficients")
    ax[0].boxplot(CS)
    ax[1].set_title("RMSD")
    ax[1].boxplot(rmsd)
    plt.savefig(i+".png")
    if "-iplot" in sys.argv:
        plt.show()


                
            

lawea=[]
for i in sys.argv[1:]:
    gmx_out_reader(i,lawea)



CompareMethods(lawea)


    
    


       
      
    
