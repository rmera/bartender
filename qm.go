/*
 * qm.go, part of Bartender
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
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
	v3 "github.com/rmera/gochem/v3"
)

//uses the ree and xtb program to run a replica-exchange simulation
func REMD(coord *v3.Matrix, mol chem.AtomMultiCharger, MD *MDSettings) string {
	if MD.method == "" {
		MD.method = "gfn0" //the default
	}
	if MD.temp == 0 {
		MD.temp = 298.0
	}
	if MD.time <= 0 {
		MD.time = 500
	}
	if MD.cpus < 0 {
		MD.cpus = runtime.NumCPU()
	}
	charge := mol.Charge()
	multi := mol.Multi()
	if MD.replicas <= 0 {
		MD.replicas = MD.cpus / 3
		if MD.replicas < 1 {
			MD.replicas = 5
		}
	}
	inpgeo := "toree.pdb"
	err := chem.PDBFileWrite(inpgeo, coord, mol, nil)
	if err != nil {
		panic(err.Error())
	}
	var reexec string = os.ExpandEnv("${BTROOT}/RE/ree")
	command := fmt.Sprintf("%s -maxt %5.3f -method %s -cpus %d -charge %d -multi %d -tref %5.3f -dielectric %3.1f -replicas %d %s %d", reexec, MD.maxtemp, MD.method, MD.cpus, charge, multi, MD.temp, MD.dielectric, MD.replicas, inpgeo, MD.time)
	LogV(2, command) /////////////////////////
	remd := exec.Command("sh", "-c", command)
	err = remd.Run()
	if err != nil {
		panic(err.Error())
	}
	return "fulltraj.xyz"

}

//This has been tested and seems to work fine.
//note that the "notused" parameter is only there to keep the same signature as REMD
func MD(coord *v3.Matrix, mol chem.AtomMultiCharger, MD *MDSettings) string {
	var dry bool = false
	if _, ok := eps2Solvent[MD.dielectric]; !ok {
		MD.dielectric = 21.0
	}
	Q := new(qm.Calc)
	if MD.method == "" {
		MD.method = "gfn0" //the default
	}
	if MD.temp == 0 {
		MD.temp = 298.0
	}
	if MD.time <= 0 {
		dry = true
		MD.time = 500 //0.5 ns, rather demanding
	}
	if MD.cpus < 0 {
		MD.cpus = runtime.NumCPU()
	}
	Q.Method = MD.method
	Q.Dielectric = MD.dielectric //acetone. Want to have an intermediate dielectric. Will be ignored for gfn0 anyway.
	Q.Job = qm.Job{MD: true}
	Q.MDTime = MD.time //simulation time (whatever unit the program uses!) it's ps for xtb
	Q.MDTemp = MD.temp
	xtb := qm.NewXTBHandle()
	xtb.SetnCPU(MD.cpus)
	err := xtb.BuildInput(coord, mol, Q)
	if err != nil {
		panic(err.Error()) //you know the drill
	}
	if !dry {
		err = xtb.Run(true) //we wait for the simulation to end, this will take a while!
	}
	//We will try to remove scoord files left by xtb, but if it doesn't work, it doesn't work
	//the program will just keep running.
	toremove, err := filepath.Glob("scoord*")
	if err == nil {
		for _, f := range toremove {
			_ = os.Remove(f)
		}
	}
	//I need to check that the thing ran correctly. There is probably some file produced a the
	//end by xtb for which I can check.
	if f, err := os.Open("xtb.trj"); err == nil {
		f.Close()        //just wanted to know if it was there. Surely there is a more direct method, but I'm lazy.
		return "xtb.trj" //the name of the resulting trajectory. It's just an multi-xyz file
	}
	panic(err.Error())
}

var eps2Solvent = map[float64]string{
	80.0: "h2o",
	5.0:  "chcl3",
	9.0:  "ch2cl2",
	21.0: "acetone",
	37.0: "acetonitrile",
	33.0: "methanol",
	2.0:  "toluene",
	7.0:  "thf",
	47.0: "dmso",
	38.0: "dmf",
	-1:   "vac",
}
