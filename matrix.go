//======================================================================
//
//        Copyright (C) 2013 Yao Peng   
//        All rights reserved
//
//        filename :matrix.go
//        description :some real matrix operations
//
//        mail:yaopengdn@126.com
//        http://weibo.com/u/2151926144
//
//======================================================================
package matrix

import (
	"fmt"
	. "math"
)

const (
	zero = 1.0e-9
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
		if Abs(tmp[maxIndex][i]) < zero {
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
		if Abs(tmp[maxIndex][i]) < zero {
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
		if Abs(tmp[maxIndex][i]) < zero {
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

func (m *Matrixf64) pivotForLU(row int) int {
	index := row
	rn := m.GetRowNum()
	sc := 0.0
	for i := row; i < rn; i++ {
		sum := 0.0
		for q := 0; q < row; q++ {
			sum += (*m)[i][q] * (*m)[q][row]
		}
		s := (*m)[i][row] - sum
		if Abs(s) > Abs(sc) {
			index = i
			sc = s
		}
	}

	if Abs(sc) < zero {
		panic("Matrix is singular to working precision.")
	}

	return index
}

func (m *Matrixf64) LU() *Matrixf64 {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	tmp := NewMatrixf64(*m)

	for i := 0; i < rn; i++ {
		tmp.pivotForLU(i)

		for j := i; j < rn; j++ {
			sum := 0.0
			for q := 0; q < i; q++ {
				sum += tmp[i][q] * tmp[q][j]
			}
			tmp[i][j] -= sum
		}

		for j := i + 1; j < rn; j++ {
			sum := 0.0
			for q := 0; q < i; q++ {
				sum += tmp[j][q] * tmp[q][i]
			}
			tmp[j][i] = (tmp[j][i] - sum) / tmp[i][i]
		}
	}
	return &tmp
}

/* Column pivoting LU decomposition */
func (m *Matrixf64) MELU() *Matrixf64 {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	tmp := NewMatrixf64(*m)

	for i := 0; i < rn; i++ {
		maxIndex := tmp.pivotForLU(i)

		if maxIndex != i {
			tmp.exchange(i, maxIndex)
		}

		for j := i; j < rn; j++ {
			sum := 0.0
			for q := 0; q < i; q++ {
				sum += tmp[i][q] * tmp[q][j]
			}
			tmp[i][j] -= sum
		}

		for j := i + 1; j < rn; j++ {
			sum := 0.0
			for q := 0; q < i; q++ {
				sum += tmp[j][q] * tmp[q][i]
			}
			tmp[j][i] = (tmp[j][i] - sum) / tmp[i][i]
		}
	}

	return &tmp
}

func (m *Matrixf64) isVector() bool {
	return len((*m)[0]) == 1
}

func (m *Matrixf64) VectorToSlice() []float64 {
	if !m.isVector() {
		panic("It must be a vector rather than a matrix.")
	}
	ret := make([]float64, len(*m))
	for i, r := range *m {
		ret[i] = r[0]
	}
	return ret
}

/*  please make sure the matrix 'm' is an LU matrix. 
for example, a non-singular equation Ax = b, the code is below:
				x := A.LU().Solve(b)								*/
func (m *Matrixf64) Solve(b *Matrixf64) *Matrixf64 {
	if !b.isVector() {
		panic("Input argument 'b' must be a vector rather than a matrix.")
	}

	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	if rn != b.GetRowNum() {
		panic("Matrix dimensions must agree.")
	}

	tmp := NewMatrixf64(*b)

	for i := 0; i < rn; i++ {
		sum := 0.0
		for j := 0; j < i; j++ {
			sum += (*m)[i][j] * (*m)[j][j] * tmp[j][0]
		}
		tmp[i][0] = (tmp[i][0] - sum) / (*m)[i][i]
	}

	for i := rn - 1; i >= 0; i-- {
		sum := 0.0
		for j := i + 1; j < rn; j++ {
			sum += (*m)[i][j] * tmp[j][0]
		}
		tmp[i][0] = tmp[i][0] - sum/(*m)[i][i]
	}
	return &tmp
}

func (m *Matrixf64) VectorNorm1() (ret float64) {
	if !m.isVector() {
		panic("It must be a vector rather than a matrix.")
	}

	for _, r := range *m {
		ret += Abs(r[0])
	}
	return
}

func (m *Matrixf64) Dist() (ret float64) {
	return m.VectorNorm2()
}

func (m *Matrixf64) VectorNorm2() (ret float64) {
	if !m.isVector() {
		panic("It must be a vector rather than a matrix.")
	}

	for _, r := range *m {
		ret += Pow(r[0], 2.0)
	}
	ret = Sqrt(ret)
	return
}

func (m *Matrixf64) VectorNormInf() (ret float64) {
	if !m.isVector() {
		panic("It must be a vector rather than a matrix.")
	}

	ret = Abs((*m)[0][0])
	for i := 1; i < m.GetRowNum(); i++ {
		if v := Abs((*m)[i][0]); v > ret {
			ret = v
		}
	}
	return
}

func (m *Matrixf64) Equal(rm *Matrixf64) bool {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != rm.GetRowNum() || cn != rm.GetColNum() {
		return false
	}

	for i := 0; i < rn; i++ {
		for j := 0; j < cn; j++ {
			if Abs((*m)[i][j]-(*rm)[i][j]) > zero {
				return false
			}
		}
	}
	return true
}

/* please make sure 'm' has enough space to accommodate the smaller matrix 'im' , or system will throw a panic*/
func (m *Matrixf64) Insert(im *Matrixf64, row, col int) *Matrixf64 {
	tmp := NewMatrixf64(*m)
	for i, r := range *im {
		for j, v := range r {
			im := i + row
			jm := j + col
			tmp[im][jm] = v
		}
	}
	return &tmp
}

/* please make sure 'm' is big enough to contain these rows or columns, or system will throw a panic*/
func (m *Matrixf64) Piece(row, col []int) *Matrixf64 {
	prn := len(row)
	pcn := len(col)

	tmp := Of64(prn, pcn)
	for _i, i := range row {
		for _j, j := range col {
			tmp[_i][_j] = (*m)[i][j]
		}
	}
	return &tmp
}

func Of64(rows, cols int) Matrixf64 {
	var m Matrixf64
	m = make([][]float64, rows)

	for i := 0; i < rows; i++ {
		m[i] = make([]float64, cols)
	}
	return m
}

func sign(d float64) int {
	switch {
	case d > zero:
		return 1
	case Abs(d) < zero:
		return 0
	case d < -zero:
		return -1
	}
	return -1
}

func (m *Matrixf64) Hess() (Q, B Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	A := NewMatrixf64(*m)
	H := Of64(rn, cn)
	I := If64(rn)
	Q = If64(rn)
	B = If64(rn)
	for r := 0; r < rn-2; r++ {
		alen := rn - r - 1
		row := make([]int, alen)

		for i := 0; i < alen; i++ {
			row[i] = r + 1 + i
		}
		col := make([]int, 1)
		col[0] = r
		a := A.Piece(row, col)
		d := a.VectorNorm2()
		var c float64
		if Abs(d) < zero {
			H = If64(rn)
		} else {
			if Abs((*a)[0][0]) < zero {
				c = d
			} else {
				c = float64(-sign((*a)[0][0])) * d
			}

			(*a)[0][0] -= c
			u := Of64(rn, 1)
			ur := u.Insert(a, r+1, 0)
			H = *I.Sub(ur.Mul(ur.T()).DotMul(2.0).DivDot((*ur.T().Mul(ur))[0][0]))
		}
		B = *H.Mul(&A).Mul(&H)
		A = B
		Q = *Q.Mul(&H)
	}
	return
}
