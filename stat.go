/*
 * stat.go, part of Bartender
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
	"math"
	"sort"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/stat"
)

//takes a slice with values (angles, distances, dihedrals) and, from their relative abundance, obtains an energy
func IBoltzmann(inp []float64, increment, temperature float64) ([]float64, []float64) {
	sort.Float64s(inp)
	divs := make([]float64, 0, 10)
	hpoints := make([]float64, 0, 10)
	for i := inp[0] - increment; i <= inp[len(inp)-1]+increment; i += increment {
		divs = append(divs, i)
		if len(divs) >= 2 {
			hpoints = append(hpoints, (divs[len(divs)-1]+divs[len(divs)-2])/2.0)
		}
	}
	//	fmt.Println(hpoints) //////////////
	histo := make([]float64, len(divs)-1)
	histo = stat.Histogram(histo, divs, inp, nil)

	//	fmt.Println(histo) /////
	//now we invert the Maxwell-Boltzmann distribution to
	energies := make([]float64, len(histo))
	largest := 0.0
	for _, v := range histo {
		if v > largest {
			largest = v
		}
	}
	for i, v := range histo {
		energies[i] = -1 * chem.R * temperature * math.Log(v/largest) //math.Log is the natural log
	}
	//we remove now the points with 0 frequency
	cleanE := make([]float64, 0, len(histo))
	cleanpoints := make([]float64, 0, len(histo))
	for i, v := range histo {
		if v != 0 {
			cleanE = append(cleanE, energies[i])
			cleanpoints = append(cleanpoints, hpoints[i])
		}
	}

	//	fmt.Println("values:", cleanpoints, "\nener:", cleanE) //////////////
	return cleanpoints, cleanE

}

type bendtor struct {
	b1  []float64
	b2  []float64
	tor []float64
}

func Newbendtor(N int) *bendtor {
	bt := new(bendtor)
	bt.tor = make([]float64, N, N)
	bt.b1 = make([]float64, N, N)
	bt.b2 = make([]float64, N, N)
	return bt
}

type bendtorFreq struct {
	tor []float64
	b1  []float64
	b2  []float64
	n   int
	e   float64
}

func NewbendtorFreq() *bendtorFreq {
	r := new(bendtorFreq)
	r.tor = make([]float64, 2, 2)
	r.b1 = make([]float64, 2, 2)
	r.b2 = make([]float64, 2, 2)
	r.n = 0 //not reaaally needed
	return r
}

type freqs []*bendtorFreq

func (f freqs) LargestFreq() int {
	var largest int
	for _, v := range f {
		if v.n > largest {
			largest = v.n
		}
	}
	return largest
}

func (f freqs) removeZeros() freqs {
	f2 := make([]*bendtorFreq, 0, len(f))
	for _, v := range f {
		if v.n != 0 {
			f2 = append(f2, v)
		}
	}
	return freqs(f2)
}

//If any of this fails, check that I didn't copy-paste the variables wront (i.e. assigned b1 to tor, for instance).

//returns the values in f but as 3 separate, matching, slices of torsions, bend1, bend2 and energies.
func (f freqs) tbbe() ([]float64, []float64, []float64, []float64) {
	t := make([]float64, 0, len(f))
	b1 := make([]float64, 0, len(f))
	b2 := make([]float64, 0, len(f))
	e := make([]float64, 0, len(f))
	for _, v := range f {
		tor := (v.tor[0] - v.tor[1]) / 2
		bend1 := (v.b2[0] - v.b1[1]) / 2
		bend2 := (v.b2[0] - v.b2[1]) / 2
		t = append(t, tor)
		b1 = append(b1, bend1)
		b2 = append(b2, bend2)
		e = append(e, v.e)
	}
	return t, b1, b2, e

}

//CountAngles will search bt and count all the elements that have the torsion and 2 angles
//within the ranges specified in AF. The final count will be loaded to the "N" field of AF
//and the same -modified- AF will be returned.
func CountAngles(AF *bendtorFreq, bt *bendtor) *bendtorFreq {
	for i, v := range bt.tor {
		if v >= AF.tor[0] && v < AF.tor[1] {
			if bt.b1[i] >= AF.b1[0] && bt.b1[i] < AF.b1[1] {
				if bt.b2[i] >= AF.b2[0] && bt.b2[i] < AF.b2[1] {
					AF.n++
				}

			}
		}

	}
	return AF
}

//the only solution here seems to be to set up a 4D slice containing the frequencies for each combination of dihedral, angle1
//and angle2. To normalize that and then transform to a 4D slice of energies, which is then used by the fitting function.

//takes a slice with values for 2 bendings and the torsion between them. From their relative abundance, obtains an energy
//the first element in increments is the increment for the torsion, the second, for the 2 angles
func IBoltzmannBT(inpt, inpb1, inpb2, incre []float64, temperature float64) ([]float64, []float64, []float64, []float64) {
	bt := Newbendtor(len(inpb1))
	copy(bt.b1, inpb1)
	copy(bt.b2, inpb2)
	copy(bt.tor, inpt)
	sort.Float64s(inpb1)
	sort.Float64s(inpb2)
	sort.Float64s(inpt)
	allFreqs := make([]*bendtorFreq, 0, 100)
	for i := inpt[0]; i <= inpt[len(inpt)-1]+incre[0]; i += incre[0] {

		for j := inpb1[0]; j <= inpb1[len(inpb1)-1]+incre[1]; j += incre[1] {

			for k := inpb2[0]; k <= inpb2[len(inpb2)-1]+incre[2]; k += incre[2] {

				freq := NewbendtorFreq()
				freq.tor = []float64{i, i + incre[0]}
				freq.b1 = []float64{j, j + incre[1]}
				freq.b2 = []float64{k, k + incre[2]}
				freq = CountAngles(freq, bt)
				allFreqs = append(allFreqs, freq)

			}
		}
	}
	//now I need to define a type for a slice of bendtorFreq with a method that gives me the total/largest "frequencies" for the set of angles.
	F := freqs(allFreqs)
	largest := F.LargestFreq()
	F = F.removeZeros()
	for _, v := range F {
		q := float64(v.n) / float64(largest)
		v.e = -1 * chem.R * temperature * math.Log(q) //math.Log is the natural log
	}
	return F.tbbe()

}
