/*
 * files.go, part of Bartender
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
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

//parses the input file, returns an slice of slices of ints,
//where the nth slice contains the atoms in the nth Martini bead
//and aslice of slices of float64, where the nth slice contains, for each atom, the fraction of the
//that atom that belonging to the nth bead. Thus, a bead can contain half an atom, for instance.
func ParseInputGeo(inpname string) (map[string][][]int, [][]int) {

	param := map[string][][]int{
		"bonds":  make([][]int, 0, 0),
		"angles": make([][]int, 0, 0),
		"reb":    make([][]int, 0, 0),
		"dihe":   make([][]int, 0, 0),
		"improp": nil,
	}
	var marked [][]int
	reading := ""
	finp, err := os.Open(inpname)
	if err != nil {
		panic(err.Error())
	}
	inp := bufio.NewReader(finp)
	for {
		line, err := inp.ReadString('\n')
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				break
			} else {
				panic(err.Error())
			}
		}
		if strings.HasPrefix(line, "#") {
			continue //comment line
		}
		if strings.HasPrefix(line, "BEADS") {
			reading = ""
			continue
		}
		if strings.HasPrefix(line, "BONDS") {
			reading = "bonds"
			continue
		}
		if strings.HasPrefix(line, "REB") {
			reading = "reb"
			continue
		}

		if strings.HasPrefix(line, "ANGLES") {
			reading = "angles"
			continue
		}
		if strings.HasPrefix(line, "DIHEDRALS") {
			reading = "dihe"
			continue
		}
		if strings.HasPrefix(line, "IMPROPERS") {
			reading = "improp"
			param["improp"] = make([][]int, 0, 0) //there may not always be impropers, so this one is only created here. Maybe I should do this for angles and dihedrals too.
			continue
		}
		if reading == "" {
			continue
		}
		pf := strings.ReplaceAll(line, " ", "")
		pf = strings.ReplaceAll(pf, "\n", "")
		if pf == "" {
			continue //shouldn't happen, but you know how users are.
		}
		star := strings.HasSuffix(pf, "*")
		pf = strings.TrimRight(pf, "*")
		fields := strings.Split(pf, ",")
		nums := make([]int, len(fields))
		for i, v := range fields {
			nums[i], err = strconv.Atoi(v)
			if err != nil {
				//		fmt.Println(fields) ///////////
				panic(err.Error())
			}
			nums[i]-- //convert from 1-based indexes to 0-based

		}
		param[reading] = append(param[reading], nums)
		if star {
			marked = append(marked, nums)
		}

	}
	fmt.Println(param) //////////////////
	finp.Close()
	return param, marked
}

//Parses the input file, returns an slice of slices of ints,
//where the nth slice contains the atoms in the nth Martini bead
//and aslice of slices of float64, where the nth slice contains, for each atom, the fraction of the
//that atom that belonging to the nth bead. Thus, a bead can contain half an atom, for instance.
func ParseInputBead(inpname string) ([][]int, [][]float64) {
	beadslice := make([][]int, 0, 0)
	wslice := make([][]float64, 0, 0)
	finp, err := os.Open(inpname)
	if err != nil {
		panic(err.Error())
	}
	inp := bufio.NewReader(finp)
	reading := false
	for {
		line, err := inp.ReadString('\n')
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				break
			} else {
				panic(err.Error())
			}
		}
		if strings.HasPrefix(line, "#") {
			continue //comment line
		}
		if strings.HasPrefix(line, "BEADS") {
			reading = true
			continue
		}
		if strings.HasPrefix(line, "BONDS") { //the "BEADS" part is over
			break
		}
		if !reading {
			continue
		}
		//now the actual reading!
		pf := strings.ReplaceAll(strings.Fields(line)[1], " ", "")
		fields := strings.Split(pf, ",")
		beadslice = append(beadslice, make([]int, len(fields)))
		wslice = append(wslice, make([]float64, len(fields)))
		m1 := len(beadslice) - 1
		for i, v := range fields {
			var bead string
			var weight string = "1.0"
			wslice[m1][i] = 1.0
			bead = v
			if strings.Contains(v, "/") {
				info := strings.Split(v, "/")
				bead = info[0]
				weight = info[1]
			}
			wslice[m1][i], err = strconv.ParseFloat(weight, 64)
			if err != nil {
				panic(err.Error())
			}
			wslice[m1][i] = 1 / wslice[m1][i]
			beadslice[m1][i], err = strconv.Atoi(bead)
			if err != nil {
				panic(err.Error())
			}
			beadslice[m1][i]-- //to convert from 1-based indexes to 0-based indexes

		}

	}
	finp.Close()
	return beadslice, wslice
}
