/*
 * mol_output.go, part of Bartender
 *
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
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/dcd"
	v3 "github.com/rmera/gochem/v3"
)

func MakePDB(coord *v3.Matrix, mol chem.Atomer, indexes [][]int) {
	binterval := 100.0 / float64(len(indexes))
	bfacs := make([]float64, mol.Len())
	for i, _ := range bfacs {
		bfacs[i] = -1
	}

	for i, v := range indexes {
		for _, w := range v {
			ifloat := float64(i)
			at := mol.Atom(w)
			at.MolID = i + 1
			if bfacs[w] == -1 {
				bfacs[w] = binterval * ifloat

			} else {
				bfacs[w] = (binterval*ifloat + bfacs[w]) / 2 //this is for atoms that are shared among beads. It will only help if the sharing beads come one after the other (say, if the atom is shared between beads 7 and 8, but not if the sharing beads are 8 and 10.
			}
		}
	}
	chem.PDBFileWrite("Beads.pdb", coord, mol, bfacs)

}

//Saves the multi-xyz trajectory in the file trajname to a DCD trajectory in the file fname.
//Since this is not really  needed, it doesn't panic. Will just return an error to be printed by main.
func DCDSave(fname, trajname string) error {
	if !strings.Contains(trajname, ".trj") {
		return fmt.Errorf("Option only available for multi-xyz trajectories, such as those produced by xtb")
	}
	mol, err := chem.XYZFileRead(trajname)
	if err != nil {
		return err
	}
	dcd, err := dcd.NewWriter(fname, mol.Len())
	coord := v3.Zeros(mol.Len())
	for {
		err = mol.Next(coord)
		if err != nil {
			break
		}
		err = dcd.WNext(coord)
		if err != nil {
			return err
		}
	}
	if _, ok := err.(chem.LastFrameError); ok { //just means we read the whole thing.
		return nil
	}
	return err

}
