package matrix

import (
	"fmt"
	. "math"
)

type Matrixf64 [][]float64

/* Matrixf64 */
func NewMatrixf64(any [][]float64) Matrixf64 {
	var m Matrixf64
	rn := len(any)
	m = make([][]float64, rn)

	for i := 0; i < rn; i++ {
		m[i] = make([]float64, len(any[i]))
		copy(m[i], any[i])
	}
	return m
}

func NewVectorf64(any []float64) Matrixf64 {
	var m Matrixf64
	rows := len(any)
	m = make([][]float64, rows)

	for i := 0; i < rows; i++ {
		m[i] = make([]float64, 1)
		m[i][0] = any[i]
	}
	return m
}

func If64(d int) Matrixf64 {
	var m Matrixf64
	m = make([][]float64, d)

	for i := 0; i < d; i++ {
		m[i] = make([]float64, d)
		m[i][i] = 1.0
	}
	return m
}

func (m *Matrixf64) Print() {
	for _, r := range *m {
		for _, c := range r {
			fmt.Printf("%.6f  ", c)
		}
		fmt.Println("")
	}
}

func (m *Matrixf64) GetVectorSum() (sum float64) {
	switch {
	case len(*m) == 1:
		for _, v := range (*m)[0] {
			sum += v
		}
	case len((*m)[0]) == 1:
		for i := 0; i < len(*m); i++ {
			sum += (*m)[i][0]

		}
	default:
		panic("It is not a vector!")
	}
	return
}

func (m *Matrixf64) T() *Matrixf64 {
	var mt Matrixf64
	mt = make([][]float64, len((*m)[0]))
	c := len(*m)

	for i, _ := range (*m)[0] {
		mt[i] = make([]float64, c)
	}

	for i, r := range *m {
		for j, v := range r {
			mt[j][i] = v
		}
	}
	return &mt
}

func (m *Matrixf64) GetRowNum() int {
	return len(*m)
}

func (m *Matrixf64) GetColNum() int {
	return len((*m)[0])
}

func (m *Matrixf64) Add(rm *Matrixf64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	mout := Matrixf64(make([][]float64, rn))
	if rn != rm.GetRowNum() || cn != rm.GetColNum() {
		panic("Matrix dimensions must agree.")
	}

	for i := 0; i < rn; i++ {
		mout[i] = make([]float64, cn)
		for j := 0; j < cn; j++ {
			mout[i][j] = (*m)[i][j] + (*rm)[i][j]
		}
	}
	return &mout
}

func (m *Matrixf64) Sub(rm *Matrixf64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	mout := Matrixf64(make([][]float64, rn))
	if rn != rm.GetRowNum() || cn != rm.GetColNum() {
		panic("Matrix dimensions must agree.")
	}

	for i := 0; i < rn; i++ {
		mout[i] = make([]float64, cn)
		for j := 0; j < cn; j++ {
			mout[i][j] = (*m)[i][j] - (*rm)[i][j]
		}
	}
	return &mout
}

func (m *Matrixf64) Mul(rm *Matrixf64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	rrn := rm.GetRowNum()
	rcn := rm.GetColNum()
	mout := Matrixf64(make([][]float64, rn))
	if cn != rrn {
		panic("Inner matrix dimensions must agree.")
	}

	for i := 0; i < rn; i++ {
		mout[i] = make([]float64, rcn)
	}

	for i := 0; i < rn; i++ {
		for k := 0; k < rcn; k++ {
			mout[i][k] = 0.0
			for j := 0; j < cn; j++ {
				mout[i][k] += (*m)[i][j] * (*rm)[j][k]
			}
		}
	}
	ret = &mout
	return ret
}

func (m *Matrixf64) exchange(i, j int) {
	for k := 0; k < m.GetColNum(); k++ {
		(*m)[i][k], (*m)[j][k] = (*m)[j][k], (*m)[i][k]
	}

}

func (m *Matrixf64) multiple(index int, mul float64) {
	for j := 0; j < m.GetColNum(); j++ {
		(*m)[index][j] *= mul
	}

}

func (m *Matrixf64) pivot(row int) int {
	index := row
	rn := m.GetRowNum()

	for i := row + 1; i < rn; i++ {
		if Abs((*m)[i][row]) > Abs((*m)[index][row]) {
			index = i
		}
	}
	return index
}

func (m *Matrixf64) multipleAdd(index, src int, mul float64) {
	for j := 0; j < m.GetColNum(); j++ {
		(*m)[index][j] += (*m)[src][j] * mul
	}
}

func (m *Matrixf64) Rank() int {
	rn := m.GetRowNum() 
	if rn > m.GetColNum() {
		return m.T().Rank()
	}
	
	tmp := NewMatrixf64(*m)
	for i := 0; i < rn; i++ {
		maxIndex := tmp.pivot(i)
		if Abs(tmp[maxIndex][i]) < 1e-9 {
			return i
		}

		if maxIndex != i {
			tmp.exchange(i, maxIndex)
		}

		for j := i + 1; j < rn; j++ {
			dMul := -tmp[j][i] / tmp[i][i]
			tmp.multipleAdd(j, i, dMul)
		}
	}
	return rn
}

func (m *Matrixf64) Det() (ret float64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	tmp := NewMatrixf64(*m)
	ret = 1
	for i := 0; i < rn; i++ {
		maxIndex := tmp.pivot(i)
		if Abs(tmp[maxIndex][i]) < 1e-9 {
			return 0.0
		}

		if maxIndex != i {
			tmp.exchange(i, maxIndex)
			ret *= -1
		}

		for j := i + 1; j < rn; j++ {
			dMul := -tmp[j][i] / tmp[i][i]
			tmp.multipleAdd(j, i, dMul)
		}
	}

	for i, _ := range tmp {
		ret *= tmp[i][i]
	}
	return ret
}

func (m *Matrixf64) Inv() *Matrixf64 {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	tmp := NewMatrixf64(*m)
	ret := If64(rn)

	for i := 0; i < rn; i++ {
		maxIndex := tmp.pivot(i)
		if Abs(tmp[maxIndex][i]) < 1e-9 {
			panic("Matrix is singular to working precision.")
		}

		if maxIndex != i {
			tmp.exchange(i, maxIndex)
			ret.exchange(i, maxIndex)
		}

		ret.multiple(i, 1.0/tmp[i][i])
		tmp.multiple(i, 1.0/tmp[i][i])

		for j := i + 1; j < rn; j++ {
			dMul := -tmp[j][i] / tmp[i][i]
			tmp.multipleAdd(j, i, dMul)
			ret.multipleAdd(j, i, dMul)
		}
	}

	for i := rn - 1; i > 0; i-- {
		for j := i - 1; j >= 0; j-- {
			dMul := -tmp[j][i] / tmp[i][i]
			tmp.multipleAdd(j, i, dMul)
			ret.multipleAdd(j, i, dMul)
		}
	}
	return &ret
}

func (m *Matrixf64) IsSymmetric() bool {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		return false
	}

	for i := 0; i < rn; i++ {
		for j := i + 1; j < cn; j++ {
			if (*m)[i][j] != (*m)[j][i] {
				return false
			}
		}
	}
	return true
}

func (m *Matrixf64) IsSquare() bool {
	return m.GetRowNum() == m.GetColNum()
}

func (m *Matrixf64) DotMul(d float64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	mout := NewMatrixf64(*m)

	for i := 0; i < rn; i++ {
		for j := 0; j < cn; j++ {
			mout[i][j] *= d
		}
	}
	return &mout
}

func (m *Matrixf64) DotDiv(d float64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	mout := NewMatrixf64(*m)

	for i := 0; i < rn; i++ {
		for j := 0; j < cn; j++ {
			mout[i][j] = d / mout[i][j]
		}
	}
	return &mout
}

func (m *Matrixf64) DivDot(d float64) (ret *Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	mout := NewMatrixf64(*m)

	for i := 0; i < rn; i++ {
		for j := 0; j < cn; j++ {
			mout[i][j] /= d
		}
	}
	return &mout
}

