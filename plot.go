package main

import (
	"image/color"
	"math"
	"strings"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

//This is mostly repeated code from fit_go. I might refactor things so these functions are called there to eliminate that problem.

func sperf(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]
	n := par[2]
	return func(x float64) float64 { return k * (1 + math.Cos(n*x-eq)) }
}

func hookef(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]
	return func(x float64) float64 { return 0.5 * k * math.Pow((x-eq), 2.0) }
}

func rebf(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]

	return func(x float64) float64 {
		//Here I also changed "eq" for "math.Cos(eq)". This seems to agree with the Gromacs
		//Manual so I think the previous behavior was a bug (see also the comment in the ReB part of the fit_go.go file)
		return 0.5 * k * math.Pow((math.Cos(x)-math.Cos(eq)), 2.0) * (1 / math.Pow(math.Sin(x), 2))
	}
}

func rybef(p []float64) func(float64) float64 {
	ret := func(x float64) float64 {
		psi := x - math.Pi
		cos := math.Cos
		pow := math.Pow
		return p[0] + p[1]*cos(psi) + p[2]*pow(cos(psi), 2) + p[3]*pow(cos(psi), 3) + p[4]*pow(cos(psi), 4) + p[5]*pow(cos(psi), 5)
	}
	return ret
}

func cosanglef(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]

	return func(x float64) float64 {
		return 0.5 * k * math.Pow((math.Cos(x)-math.Cos(eq)), 2.0)
	}
}

//plots y and f(x) vs x, as points and a line, respectively, unless given true in noplot, in which case, does nothing.
//the plot is saved to a file name.png
func Plot(f func(float64) float64, x, y []float64, name string, noplot bool) error {
	if noplot {
		return nil
	}
	name = strings.ReplaceAll(name, " ", "")
	unit_conv := 1.0
	has := strings.Contains
	if has(name, "ngle") || has(name, "Rycka") || has(name, "Cos") || has(name, "eriodic") || has(name, "mproper") {
		unit_conv = chem.Rad2Deg
	}
	pointsData := pointsPlot(x, y, unit_conv)
	funcData := funcPlot(x, f, unit_conv)
	p, err := plot.New()
	if err != nil {
		return err
	}
	p.Title.Text = name
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	// Draw a grid behind the data
	p.Add(plotter.NewGrid())
	// Make a scatter plotter and set its style.
	s, err := plotter.NewScatter(pointsData)
	if err != nil {
		return err
	}
	s.GlyphStyle.Shape = draw.PyramidGlyph{}
	s.GlyphStyle.Color = color.RGBA{R: 255, A: 255}

	// Make a line plotter and set its style.
	l, err := plotter.NewLine(funcData)
	if err != nil {
		return err
	}
	l.LineStyle.Color = color.RGBA{B: 255, A: 255}
	p.Add(s, l)
	p.Legend.Add("Trajectory", s)
	p.Legend.Add("Fitted function", l)
	// Save the plot to a PNG file.
	if err := p.Save(6*vg.Inch, 6*vg.Inch, name+".png"); err != nil {
		return err
	}

	return nil

}

func pointsPlot(x, y []float64, unit float64) plotter.XYs {
	pts := make(plotter.XYs, len(x))
	for i, v := range x {
		pts[i].X = v * unit
		pts[i].Y = y[i]
	}
	return pts
}

func funcPlot(x []float64, f func(float64) float64, unit float64) plotter.XYs {
	pts := make(plotter.XYs, len(x))
	for i, v := range x {
		pts[i].X = v * unit
		pts[i].Y = f(v)
	}
	return pts
}
