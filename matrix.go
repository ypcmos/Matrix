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
	"math/cmplx"
)

const (
	zero = 1.0e-9
)

type Matrixf64 [][]float64

type MatrixErr struct {
	err string
}

func (me *MatrixErr) Error() string {
	return me.err
}

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

func Diagf64(any []float64) Matrixf64 {
	d := len(any)
	var m Matrixf64
	m = make([][]float64, d)

	for i := 0; i < d; i++ {
		m[i] = make([]float64, d)
		m[i][i] = any[i]
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
	tmp := NewMatrixf64(*m)
	return tmp.rank()
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
func (m *Matrixf64) Replace(im *Matrixf64, row, col int) *Matrixf64 {
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

	H := Of64(rn, cn)
	I := If64(rn)
	Q = If64(rn)
	B = NewMatrixf64(*m)
	for r := 0; r < rn-2; r++ {
		alen := rn - r - 1
		row := make([]int, alen)

		for i := 0; i < alen; i++ {
			row[i] = r + 1 + i
		}
		col := make([]int, 1)
		col[0] = r
		a := B.Piece(row, col)
		d := a.VectorNorm2()
		var c float64
		ar := (*a)[0][0]
		(*a)[0][0] = 0
		if a.VectorNormInf() < zero {
			H = If64(rn)
		} else {
			if Abs(ar) < zero {
				c = d
			} else {
				c = float64(-sign(ar)) * d
			}

			(*a)[0][0] = ar - c
			u := Of64(rn, 1)
			ur := u.Replace(a, r+1, 0)
			H = *I.Sub(ur.Mul(ur.T()).DotMul(2.0).DivDot((*ur.T().Mul(ur))[0][0]))
		}
		B = *H.Mul(&B).Mul(&H)
		Q = *Q.Mul(&H)
	}
	return
}

func (m *Matrixf64) QR() (Q, R Matrixf64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	R = NewMatrixf64(*m)
	H := Of64(rn, cn)
	I := If64(rn)
	Q = If64(rn)
	for r := 0; r < rn-1; r++ {
		alen := rn - r
		row := make([]int, alen)

		for i := 0; i < alen; i++ {
			row[i] = r + i
		}
		col := make([]int, 1)
		col[0] = r
		a := R.Piece(row, col)
		d := a.VectorNorm2()
		ar := (*a)[0][0]
		var c float64
		if a.VectorNormInf() < zero {
			H = If64(rn)
		} else {
			if Abs(ar) < zero {
				c = d
			} else {
				c = float64(-sign(ar)) * d
			}

			(*a)[0][0] = ar - c
			u := Of64(rn, 1)
			ur := u.Replace(a, r, 0)
			H = *I.Sub(ur.Mul(ur.T()).DotMul(2.0).DivDot((*ur.T().Mul(ur))[0][0]))
		}
		R = *H.Mul(&R)
		Q = *Q.Mul(&H)
	}

	/* 	
		Makes 'R' with diagonal elements which are non-negative real numbers 
		A = qr = qI'Ir = qIIr = QR  =>  Q = qI, R=Ir
		(I = [...
				-1
					...
						1
							...])
	*/
	for i, _ := range R {
		if sign(R[i][i]) == -1 {
			I[i][i] = -1
		}
	}
	R = *I.Mul(&R)
	Q = *Q.Mul(&I)
	return
}

func (m *Matrixf64) eig2() (ret []complex128) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	switch rn {
	case 1:
		ret = make([]complex128, 1)
		ret[0] = complex((*m)[0][0], 0)
	case 2:
		ret = make([]complex128, 2)
		ret[0] = 0.5 * (complex((*m)[0][0]+(*m)[1][1], 0) +
			cmplx.Sqrt(complex(Pow((*m)[0][0]+(*m)[1][1], 2.0)-4*((*m)[0][0]*(*m)[1][1]-(*m)[0][1]*(*m)[1][0]), 0)))
		ret[1] = 0.5 * (complex((*m)[0][0]+(*m)[1][1], 0) -
			cmplx.Sqrt(complex(Pow((*m)[0][0]+(*m)[1][1], 2.0)-4*((*m)[0][0]*(*m)[1][1]-(*m)[0][1]*(*m)[1][0]), 0)))
	default:
		panic("Only support matirx that dimension <= 2.")
	}
	return
}

/* WARNING: this function will damage the matrix 'm' */
func (m *Matrixf64) eig(ret *[]complex128) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}
	getChunk := func(m *Matrixf64, i int) (rm *Matrixf64) {
		size := rn - i
		rc := make([]int, size)

		for _i := 0; _i < size; _i++ {
			rc[_i] = _i + i
		}
		rm = m.Piece(rc, rc)

		rc = make([]int, i)
		for _i := 0; _i < i; _i++ {
			rc[_i] = _i
		}
		*m = *m.Piece(rc, rc)
		return
	}

	qr2 := func(m *Matrixf64) Matrixf64 {
		//rd >= 3
		rd := m.GetRowNum()
		del := (*m)[rd-1][rd-1] + (*m)[rd-2][rd-2]
		rou := (*m)[rd-1][rd-1]*(*m)[rd-2][rd-2] - (*m)[rd-1][rd-2]*(*m)[rd-2][rd-1]
		p := (*m)[0][0]*((*m)[0][0]-del) + (*m)[0][1]*(*m)[1][0] + rou
		q := (*m)[1][0] * ((*m)[0][0] + (*m)[1][1] - del)
		r := (*m)[1][0] * (*m)[2][1]
		getP := func(i int, p, q, r float64) (ret Matrixf64) {
			a := NewVectorf64([]float64{p, q, r})
			if rd-i < 4 {
				a = a[0:2]
			}
			var c float64
			if Abs(p) < zero {
				c = a.VectorNorm2()
			} else {
				c = float64(sign(p)) * a.VectorNorm2()
			}
			v := Of64(rd, 1)
			a[0][0] += c

			v = *v.Replace(&a, i+1, 0)
			ret = If64(rd)
			if v.VectorNormInf() > zero {
				ret = *ret.Sub(v.Mul(v.T()).DotMul(2 / Pow(v.VectorNorm2(), 2.0)))
			}
			return
		}
		P := getP(-1, p, q, r)
		*m = *P.Mul(m).Mul(&P)
		for i := 0; i < rd-2; i++ {
			p = (*m)[i+1][i]
			q = (*m)[i+2][i]
			if i > rd-4 {
				r = 0
			} else {
				r = (*m)[i+3][i]
			}
			P = getP(i, p, q, r)
			*m = *P.Mul(m).Mul(&P)
		}
		return *m
	}
	if rn <= 2 {
		*ret = append(*ret, m.eig2()...)
	} else {
		for i := rn - 1; i > 0; i-- {
			if Abs((*m)[i][i-1]) < zero {
				sd := rn - i
				sm := getChunk(m, i)
				switch {
				case sd <= 2:
					sm.eig(ret)
				default:
					*sm = qr2(sm)
					sm.eig(ret)

				}
				i = m.GetRowNum() - 1
			}
		}
		switch {
		case m.GetRowNum() <= 2:
			m.eig(ret)
		default:
			*m = qr2(m)
			m.eig(ret)
		}
	}
	return
}

func (m *Matrixf64) Eig() (ret []complex128) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	_, tmp := m.Hess()
	tmp.eig(&ret)
	return
}

func (m *Matrixf64) Gauss() *Matrixf64 {
	tmp := NewMatrixf64(*m)
	tmp.gauss()
	return &tmp
}

/* WARNING: this function will damage the matrix 'm' */
func (m *Matrixf64) gauss() {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	d := func(i, j int) (ret int) {
		if i >= j {
			ret = j
		} else {
			ret = i
		}
		return ret
	}(rn, cn)

	getChunk := func(m *Matrixf64, i, j int) (sm *Matrixf64) {
		rows := rn - i
		cols := cn - j
		rindex := make([]int, rows)
		cindex := make([]int, cols)
		for _i := 0; _i < rows; _i++ {
			rindex[_i] = _i + i
		}

		for _i := 0; _i < cols; _i++ {
			cindex[_i] = _i + j
		}
		sm = m.Piece(rindex, cindex)
		return
	}

	for i := 0; i < d; i++ {
		maxIndex := m.pivot(i)
		if maxIndex != i {
			m.exchange(i, maxIndex)
		}
		if Abs((*m)[i][i]) < zero {
			below := getChunk(m, i, i)
			if rm := Of64(below.GetRowNum(), below.GetColNum()); below.Equal(&rm) {
				break
			}
			sm := getChunk(m, i, i+1)
			sm.gauss()
			*m = func(tm *Matrixf64) Matrixf64 {
				return *m.Replace(tm, i, i+1)
			}(sm)
			break
		}

		m.multiple(i, 1.0/(*m)[i][i])
		for j := i + 1; j < d; j++ {
			dMul := -(*m)[j][i] / (*m)[i][i]
			m.multipleAdd(j, i, dMul)

		}
	}

	if rn > d {
		x := rn - d
		O := Of64(x, cn)
		*m = *m.Replace(&O, d, 0)
	}
}

/* WARNING: this function will damage the matrix 'm' */
func (m *Matrixf64) rank() int {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	rank := 0
	if rn > cn {
		return m.T().rank()
	}

	getChunk := func(m *Matrixf64, i, j int) (sm *Matrixf64) {
		rows := rn - i
		cols := cn - j
		rindex := make([]int, rows)
		cindex := make([]int, cols)
		for _i := 0; _i < rows; _i++ {
			rindex[_i] = _i + i
		}

		for _i := 0; _i < cols; _i++ {
			cindex[_i] = _i + j
		}
		sm = m.Piece(rindex, cindex)
		return
	}

	tmp := *m
	for i := 0; i < rn; i++ {
		maxIndex := tmp.pivot(i)
		if maxIndex != i {
			tmp.exchange(i, maxIndex)
		}
		if Abs(tmp[i][i]) < zero {
			below := getChunk(&tmp, i, i)
			if rm := Of64(below.GetRowNum(), below.GetColNum()); below.Equal(&rm) {
				break
			}
			sm := getChunk(&tmp, i, i+1)
			rank += sm.rank()
			break
		} else {
			rank++
		}

		for j := i + 1; j < rn; j++ {
			dMul := -tmp[j][i] / tmp[i][i]
			tmp.multipleAdd(j, i, dMul)
		}
	}
	return rank
}

func (m *Matrixf64) MatrixNorm1() (ret float64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}
	ret = -1
	for j := 0; j < cn; j++ {
		sum := 0.0
		for i := 0; i < rn; i++ {
			sum += (*m)[i][j]
		}

		if ret < sum {
			ret = sum
		}
	}
	return
}

func (m *Matrixf64) Rou() (ret float64) {
	lenbda := m.Eig()
	ret = -1

	for _, v := range lenbda {
		abs := cmplx.Abs(v)

		if ret < abs {
			ret = abs
		}
	}
	return
}

func (m *Matrixf64) MatrixNorm2() (ret float64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	return Sqrt(m.T().Mul(m).Rou())
}

func (m *Matrixf64) MatrixNormInf() (ret float64) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}
	ret = -1
	for _, r := range *m {
		sum := 0.0
		for _, v := range r {
			sum += Abs(v)
		}

		if ret < sum {
			ret = sum
		}
	}
	return
}

func (m *Matrixf64) Cond() (ret float64) {
	return m.Inv().MatrixNorm2() * m.MatrixNorm2()
}

/*
	If you can guarantee that the coefficient matrix is nonsingular, 
	should use 'Solve' instead of this for its low efficiency
*/
func (m *Matrixf64) SolveEquation(b *Matrixf64) (x *Matrixf64, e error) {
	rn := m.GetRowNum()
	cn := m.GetColNum()
	if rn != cn {
		panic("Matrix must be square.")
	}

	tmp := m.InsertCols(b, m.GetColNum())
	if tmp.Rank() > m.Rank() {
		e = &MatrixErr{"Equation has no solution"}
		x = nil
		return
	}
	e = nil
	gauss := *tmp.Gauss()
	d := rn
	findAndExchange := func(i int) bool {
		for j := i + 1; j < d; j++ {
			if Abs(gauss[i][j]) > zero {
				gauss[i], gauss[j] = gauss[j], gauss[i]
				return true
			}
		}
		return false
	}

	for i := 0; i < d; i++ {
		if Abs(gauss[i][i]) < zero {
		again:
			if findAndExchange(i) {
				goto again
			}
			gauss[i][i] = 1
			gauss[i][d] = 1
		}
	}
	x = func() *Matrixf64 {
		row := make([]int, d)
		for i := 0; i < d; i++ {
			row[i] = i
		}
		col := []int{d}
		tmp := *gauss.Piece(row, col)

		for i := d - 1; i >= 0; i-- {
			sum := 0.0
			for j := i + 1; j < d; j++ {
				sum += gauss[i][j] * tmp[j][0]
			}
			tmp[i][0] = tmp[i][0] - sum/gauss[i][i]
		}

		sum := tmp.VectorNorm2()
		tmp = *tmp.DivDot(sum)
		return &tmp
	}()
	return
}

func (m *Matrixf64) InsertCols(im *Matrixf64, col int) *Matrixf64 {
	mrn := m.GetRowNum()
	mcn := m.GetColNum()
	irn := im.GetRowNum()
	icn := im.GetColNum()

	if mrn != irn {
		panic("The two matrix must same numbers of rows.")
	}

	rn := mrn
	cn := mcn + icn
	tmp := Of64(rn, cn)
	tmp = *tmp.Replace(im, 0, col)
	for i := 0; i < mrn; i++ {
		for j := 0; j < mcn; j++ {
			var _j int
			if j >= col {
				_j += j + icn
			} else {
				_j = j
			}
			tmp[i][_j] = (*m)[i][j]

		}
	}
	return &tmp
}

func (m *Matrixf64) InsertRows(im *Matrixf64, row int) *Matrixf64 {
	mrn := m.GetRowNum()
	mcn := m.GetColNum()
	irn := im.GetRowNum()
	icn := im.GetColNum()

	if mcn != icn {
		panic("The two matrix must same numbers of columns.")
	}

	rn := mrn + irn
	cn := mcn
	tmp := Of64(rn, cn)
	tmp = *tmp.Replace(im, row, 0)
	for i := 0; i < mrn; i++ {
		for j := 0; j < mcn; j++ {
			var _i int
			if i >= row {
				_i += i + irn
			} else {
				_i = i
			}
			tmp[_i][j] = (*m)[i][j]

		}
	}
	return &tmp
}
