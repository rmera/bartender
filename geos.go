/*
 * geos.go, part of Bartender
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
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

//Distance *in nm*!!!!!!
func distance(a, b, aux *v3.Matrix) float64 {
	aux.Sub(a, b)
	A2nm := 1 / 10.0
	return aux.Norm(2) * A2nm
}

func TrajAn(traj chem.Traj, mol chem.Atomer, indexes [][]int, weights [][]float64, wanted map[string][][]int) map[string][][]float64 {
	ret := map[string][][]float64{
		"bonds":  make([][]float64, 0),
		"angles": make([][]float64, 0),
		"reb":    make([][]float64, 0),
		"dihe":   make([][]float64, 0),
		"improp": nil,
	}
	var err error
	coord := v3.Zeros(traj.Len())
	for {
		err = traj.Next(coord)
		if err != nil {
			break
		}
		tmap := FramePar(coord, mol, indexes, weights, wanted)
		ret = UpdateMap(tmap, ret)
	}
	if _, ok := err.(chem.LastFrameError); ok { //just means we read the whole thing.
		return ret
	}
	panic(err.Error()) //something wrong happened, as the error is not LastFrameError

}

func UpdateMap(temp map[string][]float64, mmap map[string][][]float64) map[string][][]float64 {
	/*
		in temp we have, for one frame, temp["angle"]{a1,a2,a3,a4}
		in mmap we have mmap["angle"]{{a11,a12,a13...}{a21,a22,a23...}}
		so I need to append, say temp["angle"][1] to mmap[angle][1]
	*/
	for k, v := range temp {
		if v == nil || len(v) == 0 {
			continue
		}
		for l, w := range v {
			if len(mmap[k]) == 0 {
				for _, _ = range v {
					mmap[k] = append(mmap[k], make([]float64, 0, 100)) //this should happen if we are on the first frame. Note that len(mmap[k]) can never be smaller than l-2
				}
			}
			//	fmt.Println(temp, mmap, k, l) ////////////////////////////
			mmap[k][l] = append(mmap[k][l], w)
		}
	}
	return mmap
}

//Builds a map with al the wanted data in one frame
func FramePar(coord *v3.Matrix, mol chem.Atomer, indexes [][]int, weights [][]float64, wanted map[string][][]int) map[string][]float64 {
	beads := make([]*v3.Matrix, len(indexes))
	for i, v := range indexes {
		beads[i] = WCOM(coord, mol, v, weights[i])
	}
	aux := v3.Zeros(1)
	aux2 := v3.Zeros(1)

	ret := map[string][]float64{
		"bonds":  make([]float64, 0),
		"angles": make([]float64, 0),
		"reb":    make([]float64, 0),
		"dihe":   make([]float64, 0),
		"improp": nil,
	}
	//Distances first
	for _, v := range wanted["bonds"] {
		//	LogV(true, ret["bonds"])                                                     // distance(beads[v[0]], beads[v[1]], aux), v[0], v[1])              //////////////////////////////////////////////////////////////////////////
		ret["bonds"] = append(ret["bonds"], distance(beads[v[0]], beads[v[1]], aux)) //we don't check that v has the correct lenght. You are on your own there.
	}
	for _, v := range wanted["angles"] {
		aux.Sub(beads[v[0]], beads[v[1]])
		aux2.Sub(beads[v[2]], beads[v[1]])
		ret["angles"] = append(ret["angles"], chem.Angle(aux, aux2))
	}
	for _, v := range wanted["reb"] {
		aux.Sub(beads[v[0]], beads[v[1]])
		aux2.Sub(beads[v[2]], beads[v[1]])
		ret["reb"] = append(ret["reb"], chem.Angle(aux, aux2))
	}

	for _, v := range wanted["dihe"] {
		ret["dihe"] = append(ret["dihe"], chem.DihedralAlt(beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]]))
	}
	if wanted["improp"] == nil {
		return ret
	}
	ret["improp"] = make([]float64, 0)
	for _, v := range wanted["improp"] {
		ret["improp"] = append(ret["improp"], chem.Improper(beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]]))
	}
	return ret
}

//Obtains the COM for a subset of atoms from mol/coord given by indexes, where the massess are additionally
//weighted by the weights slice, if given.
//this is not a very optimized function, it allocates a lot.
func WCOM(coord *v3.Matrix, mol chem.Atomer, indexes []int, weights []float64) *v3.Matrix {
	ncoord := v3.Zeros(len(indexes))
	nmol := chem.NewTopology(0, 1)
	ncoord.SomeVecs(coord, indexes)
	nmol.SomeAtoms(mol, indexes)
	//	masses := nmol.Len()
	//	if err != nil {
	//		panic(err.Error())
	//	}
	//hopefully we always get this matrix, so we don't have to allocate a new every time.
	if weights == nil {
		weights = make([]float64, nmol.Len())
		for i, _ := range weights {
			weights[i] = 1
		}
	}
	//	massD := mat.NewDense(1, len(masses), masses)
	wD := mat.NewDense(len(weights), 1, weights)
	//	fmt.Println(len(weights), ncoord.NVecs()) ////////////////////////////////////////////////////
	//	massD.MulElem(massD, wD)
	massD := wD // Martini recommends that the centroid is used, not the COM

	gr, _ := ncoord.Dims()
	//	fmt.Println(gr) //////////////////////////////////////////////////////////////////////////////////
	tmp2 := make([]float64, gr, gr)
	for i, _ := range tmp2 {
		tmp2[i] = 1
	}
	gnOnesvector := mat.NewDense(1, gr, tmp2) //gnOnes(1, gr)
	ref := v3.Zeros(gr)
	ref.ScaleByCol(ncoord, massD)
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mat.Sum(massD), ref2)
	return ref2
}
