# This is Bartender, beta.


## Install 

To install the binary distribution, uncompress the `tgz` file
to some suitable directory, set the variable BTROOT to that
directory, and put the bartender excecutable in the PATH.

To use Bartender, you will also need the [xtb program](https://github.com/grimme-lab/xtb) from the Grimme group.

## Bartender use

Assuming the Bartender and xtb excecutables are in the PATH, you
can use Bartender by:

```
bartender [flags] Geometry.xyz BartenderInput.inp
```

The geometry file can also be in PDB or GRO format. The Bartender input format
is a simple way of specifying the bonded parameters to be obtained. You can find a
sample in the directory of the distribution, under _samples/_.

The optional flags control the way Bartender behaves. Sensible defaults have been prepared so, in
most cases, no flags are needed. Use

```
bartender -help
```

To get all the flags available and their use. We document here some of the most common flags:

*  `-charge` _int_ the total charge of the system, in a.u. (default 0)
*  `-method` _string_ The method employed in the semiempirical simulation. Valid options are gfn0, gfn1,gfn2 and gfnff (default "gfnff")
*  `-time` _int_ The total simulation time for the QM MD, in ps. If a number <0 is given, the MD will not be performed, and a previous trajectory will be used (default 1000).
*  `-verbose` _int_  Sets the level of verbosity
*  `-dcdSave` _filename.dcd_ Saves the trajectory produced by xtb in the more compact DCD format
*  `-owntraj` _filename_ Reads a trajectory (DCD, XTC, multiPDB or multiXYZ, identified from the file extension) instead of performing an xtb simulation.


## Latest changes:

1. A REMD can be requested by placing a star (`*`) near any dihedral in the
Bartender input file.
2. New flags to control the REMD
3. Instead of "verbose" or "non-verbose" there are now 4 levels of verbosity, 
from 0 to 3. The default is 1 (on the quiet side)
4. The default method is now gfnff (it could still change)


## REMD

The REMD works by calling another program, the Replica Exchange Engine (REE)
which will is open source and on [Github](https://github.com/rmera/ree), 
but is also distributed together with bartender as a binary (in the 
sub-directory RE). You can try the binary on its own (without Bartender)
to run REMDs on a given geometry (xtb is requited). It is used as:

```
./ree OPTIONS geometry.xyz SIMULATION_TIME
```

Use:

```
./ree -help
```
To see the avaliable options



