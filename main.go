/*
 * main.go, part of Bartender
 *
 *
 * Copyright 2020 Raul Mera <rmera{at}usach(dot)cl>
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 */

/*To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche*/

package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/stat"
)

//Global variables... Sometimes, you gotta use'em
var verb int

//If level is larger or equal, prints the d arguments to stderr
//otherwise, does nothing.
func LogV(level int, d ...interface{}) {
	if level <= verb {
		fmt.Fprintln(os.Stderr, d...)
	}

}

func main() {
	//There will be _tons_ of flags, but they are meant not to be needed the 99% of the time.
	avsasaskip := flag.Int("avsasaskip", 1, "If averaged-SASAs are requested, read only every nth frames of the trajectory")
	cpus := flag.Int("cpus", -1, "the total CPUs used for the QM calculations. If a number <0 is given, all logical CPUs are used")
	refit := flag.Bool("refit", false, "Only do a re-fit for the bonded parameters from an existing trajectory. Equivalent to -time 1 -nobeadtype")
	noplot := flag.Bool("noplot", false, "Do not produce the plots that would normally be written for each parameter fitted")
	owntraj := flag.String("owntraj", "", "Use the given trajectory for geometry analysis, instead of obtaining a GFN one. DCD, multi-PDB and multi-XYZ formats are allowed. XTC is allowed if the xdrfile library is installed")
	verbose := flag.Int("verbose", 0, "Print lots of additional information (mostly for debugging)")
	mdtime := flag.Int("time", 1000, "the total simulation time, in ps. If a number <0 is given, the MD will not be performed, and a previous trajectory will be used (the program will crash if no such previous trajectory is present)")
	charge := flag.Int("charge", 0, "the total charge of the system, in a.u. Needed for partial charges calculation")
	dcdsave := flag.String("dcdSave", "", "If given, Bartender will save the xtb-calculated trajectory in DCD format with the filename given")
	method := flag.String("method", "gfnff", "The method employed in the semiempirical simulation. Valid options are gfn0, gfn1,gfn2 and gfnff")
	temperature := flag.Float64("temperature", 298, "The temperature for the MD simulation, in K")
	bi := flag.Float64("bondIncrement", 0.001, "The bin width for the bond distance histograms, in nm")
	ai := flag.Float64("angleIcrcement", 1, "The bin width for the A-B-C angle histograms, in degrees")
	di := flag.Float64("dihedralIncrement", 10, "The bin width for the A-B-C-D dihedral histograms, in degrees")
	ii := flag.Float64("improperIncrement", 1, "The bin width for the A-B-C-D improper dihedral histograms, in degrees")
	dielectric := flag.Float64("dielectric", 21.0, "The dielectric constant for continuum solvent in QM calculations. Only some values are allowed (see code) -1 for vacuum calculations. Default is acetone")
	replicas := flag.Int("replicas", 0, "Number of replicas in a replica-exchange MD simulation, if performed. If less or equal zero, Bartender will come up with a reasonable number")
	maxtemp := flag.Float64("maxtemp", 400.0, "The maximum temperature allowed for a replica in a replica-exchange simulation, if performed")
	exfreq := flag.Int("exfreq", 2, "The frequency of attempted replica exchanges in a replica-exchange simulation, if performed")

	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s: [flags] geomtry.pdb/.gro/.xyz bartender_input.inp \n\nFlags:\n", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()
	verb = *verbose
	if *refit {
		*mdtime = -1
	}
	//just in case.
	if *avsasaskip < 1 {
		*avsasaskip = 1
	}
	//the angle increments will be in radians
	d2r := chem.Deg2Rad
	increments := map[string]float64{
		"bonds":  *bi,
		"angles": *ai * d2r,
		"reb":    *ai * d2r,
		"dihe":   *di * d2r,
		"improp": *ii * d2r,
	}

	param := map[string][]*bonded{
		"bonds":  make([]*bonded, 0, 0),
		"angles": make([]*bonded, 0, 0),
		"reb":    make([]*bonded, 0, 0),
		"dihe":   make([]*bonded, 0, 0),
		"improp": nil,
	}
	args := flag.Args()
	geoname := args[0]
	inpname := args[1]
	fmt.Printf("Use:\n  $BARTENDERPATH/bartender  [FLAGS] geometry_file input_file\n Use \"bartender -help\" to see the available flags\n")
	var mol *chem.Molecule
	var err error
	extension := strings.ToLower(strings.Split(geoname, ".")[1])
	switch extension {
	case "gro":
		mol, err = chem.GroFileRead(geoname)
	case "pdb":
		mol, err = chem.PDBFileRead(geoname, false)
	default:
		mol, err = chem.XYZFileRead(geoname)
	}
	if err != nil {
		panic(err.Error())
	}
	mol.SetCharge(*charge) //needed for the MD and the partial charges calculation
	wanted, marked := ParseInputGeo(inpname)
	MDEngine := MD
	if len(marked) != 0 {
		MDEngine = REMD
	}
	beads, weights := ParseInputBead(inpname)
	LogV(2, wanted, "beads:", beads, "weights:", weights)
	//bonded parameters
	MakePDB(mol.Coords[0], mol, beads)
	MDS := &MDSettings{time: *mdtime, method: *method, temp: *temperature, dielectric: *dielectric, cpus: *cpus, replicas: *replicas, maxtemp: *maxtemp, exfreq: *exfreq}

	//Here we run the calculation to get a GFN0/2 trajectory, or we read whatever trajectory the user wants to supply
	var mdout chem.Traj
	var datamap map[string][][]float64
	if true { // This "if" is for functionality that we removed temporarily, so I prefer to keep it there. Sorry about that :-)
		if *owntraj == "" {
			mdoutname := MDEngine(mol.Coords[0], mol, MDS) //This will take a while
			_, mdout, err = chem.XYZFileAsTraj(mdoutname)
			if err != nil {
				panic(err.Error())
			}
		} else {
			mdout, err = OpenTraj(*owntraj)
			if err != nil {
				panic(err.Error())
			}
		}
		datamap = TrajAn(mdout, mol, beads, weights, wanted)
	}
	//select which functions will be used, Go or Python
	HookeFit := GoHookeFit
	SimplePeriodicFit := GoSimplePeriodicFit
	RyckBelleFit := GoRyckBelleFit
	CosAngleFit := GoCosAngleFit
	ReBFit := GoReBFit
	fmt.Printf("All Energies in kJ/mol, distances in nm, angles in degrees\n")
	for k, v := range datamap {
		for i, w := range v {
			mean := stat.Mean(w, nil)
			LogV(2, "\n***", k, wanted[k][i], mean, k, len(w), increments[k])
			points, E := IBoltzmann(w, increments[k], *temperature) //doesn't return anything for now, but prints intermediate data.
			LogV(3, "Points:", points, "\nEnergies:", E)
			category := CategoryName(k)
			beadst := BeadsText(wanted[k][i])
			switch k {
			case "dihe":
				par, R2 := SimplePeriodicFit(points, E)
				if par[0] < 0 {
					LogV(1, "Dihedral corrected. Was", par[0], "became", 2*math.Pi+par[0], "Periodicity was", par[2], "Absolute value has been taken")
					par[0] = 2*math.Pi + par[0]
					par[2] = math.Abs(par[2])
				}

				LogV(3, Plot(sperf(par), points, E, fmt.Sprintf("Simple_periodic_%s", beadst), *noplot))
				par[0] = par[0] * chem.Rad2Deg
				comment := false
				if R2 > 10 {
					comment = true
				}
				b := NewBonded(i, wanted[k][i], par, R2, 1, comment)
				param[k] = append(param[k], b)
				LogV(1, fmt.Sprintf("S. Periodic. fit for the %s between  beads %s: eq: %5.3f k: %5.3f n: %5.3f Fit RMSD: %5.3f\n", category, beadst, par[0], par[1], par[2], R2))
				par2, R22 := RyckBelleFit(points, E)
				LogV(2, Plot(rybef(par2), points, E, fmt.Sprintf("Ryckaert-Belleman_%s", beadst), *noplot))
				b = NewBonded(i, wanted[k][i], par2, R22, 3, !comment) //so if simple periodic was commented, whis will not, and viceversa.

				LogV(1, fmt.Sprintf("Ryckaert-Bellemans fit for the %s between  beads %s: C1: %5.3f C2: %5.3f C3: %5.3f C4 %3.5f C5 %3.5f Fit RMSD: %5.3f\n", category, beadst, par2[0], par2[1], par2[2], par2[3], par2[4], R22))
				param[k] = append(param[k], b) //[len(param[k])-1] = append(param[k][len(param[k])-1], par2...) //just one after the other
				ia := increments["angles"]

				par3, R23 := ManageBendingTorsion(datamap, wanted, i, *temperature, []float64{increments["dihe"], ia, ia})
				//R23 should never be negative, so we'll use a negative value to signal that the fit was not obtained.
				if R23 >= 0 {
					//Ill add something to the log later -_-
					b = NewBonded(i, wanted[k][i], par3, R23, 11, true)
					param[k] = append(param[k], b)
				} else {
					LogV(1, fmt.Sprintf("Combined bending-torsion potential for beands %s will not be obtained, for lack of bending angles in input", beadst))
				}
			case "improp":
				par, R2 := HookeFit(points, E)
				Plot(hookef(par), points, E, fmt.Sprintf("Improper_Hooke_%s", beadst), *noplot)
				//par = append(par, R2)
				LogV(1, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadst, par[0]*chem.Rad2Deg, par[1], R2))
				par[0] = par[0] * chem.Rad2Deg
				impropTolerance := 10.0 //degrees
				if math.Abs(par[0]-180) < impropTolerance || math.Abs(par[0]) < impropTolerance {
					LogV(1, "**The previous equilibrium angle will be set to 180 deg! Check that it is close enough to that value, or to 0\n")
					par[0] = 180 //we force the improper dihedrals to be 180 degrees, to avoid a discontinuity in some of the functions
				} else {
					LogV(1, "**The previous equilibrium angle is too far from 180 or 0 to set it to 180 so it will")
					LogV(1, "be left as-is. This could cause numerical problems in some functions. Check that it is what you want\n")
				}
				b := NewBonded(i, wanted[k][i], par, R2, 2, false)

				param[k] = append(param[k], b)

			case "angles":
				par, R2 := HookeFit(points, E)
				LogV(3, Plot(hookef(par), points, E, fmt.Sprintf("Angle_Hooke_%s", beadst), *noplot))
				par = append(par, R2)

				par[0] = par[0] * chem.Rad2Deg
				b := NewBonded(i, wanted[k][i], par, R2, 1, false)
				param[k] = append(param[k], b)
				LogV(1, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadst, par[0], par[1], R2))

				par2, R22 := CosAngleFit(points, E)
				par2[0] = par2[0] * chem.Rad2Deg
				LogV(2, Plot(cosanglef(par2), points, E, fmt.Sprintf("CosAngle_%s", beadst), *noplot))
				b = NewBonded(i, wanted[k][i], par2, R22, 2, true)
				LogV(1, fmt.Sprintf("Cosine Angle (Gromos96) fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadst, par2[0], par2[1], R22))
				param[k] = append(param[k], b) //,[len(param[k])-1] = append(param[k][len(param[k])-1], par2...) //just one after the other

			case "bonds":
				par, R2 := HookeFit(points, E)
				LogV(3, Plot(hookef(par), points, E, fmt.Sprintf("Bond_Hooke_%s", beadst), *noplot))
				b := NewBonded(i, wanted[k][i], par, R2, 1, false)
				param[k] = append(param[k], b)
				LogV(1, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadst, par[0], par[1], R2))
			case "reb":
				par, R2 := ReBFit(points, E)
				LogV(3, Plot(rebf(par), points, E, fmt.Sprintf("ReB_%s", beadst), *noplot))
				par = append(par, R2)
				par[0] = par[0] * chem.Rad2Deg
				b := NewBonded(i, wanted[k][i], par, R2, 10, false)
				param[k] = append(param[k], b)
				LogV(1, fmt.Sprintf("Reb fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadst, par[0], par[1], R2))
			}
		}
	}
	PrintBonded(param, "gmx_out.itp")
	//	fmt.Println(datamap) //this needs to be replaced by whatever statistical analysis used to extract eq values/force constants from the distributions  in datamap

	if *owntraj == "" && *dcdsave != "" {
		err = DCDSave(*dcdsave, "xtb.trj")
		if err != nil {
			LogV(0, "Couldn't save DCD trajectory: ", err.Error()) //This always gets printed, even in non-verbose mode.
		}
	}

	fmt.Println("\nYour Martini, Mr. Bond.")

}

func BeadsText(beads []int) string {
	ret := " "
	for _, v := range beads {
		ret = ret + strconv.Itoa(v+1) + "-"
	}
	return ret[:len(ret)-1]
}

func CategoryName(k string) string {
	category := strings.ReplaceAll(k, "s", "")
	//	category = strings.Title(category)
	if category == "dihe" {
		category = "dihedral"
	} else if category == "improp" {
		category = "improper dihedral"
	}
	return category

}

type bonded struct {
	ID        int
	beads     []int
	params    []float64
	rmsd      float64
	functype  int
	commented bool
}

func NewBonded(ID int, beads []int, params []float64, rmsd float64, functype int, commented bool) *bonded {
	ret := new(bonded)
	ret.beads = beads
	ret.params = params
	ret.rmsd = rmsd
	ret.functype = functype
	ret.commented = commented
	ret.ID = ID
	return ret

}

func (b *bonded) Comment() string {
	if b.commented {
		return ";;"
	}
	return ""
}

//Settings for MD. Not all these are
//always needed.
type MDSettings struct {
	time       int
	method     string
	temp       float64
	dielectric float64
	cpus       int
	replicas   int
	maxtemp    float64
	exfreq     int
}
