package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

const (
	latX   = 1.291   // 系のx軸長さ
	latY   = 1.291   // 系のy軸長さ
	lx     = 5       // x方向格子点数
	ly     = 5       // y方向格子点数
	natom  = lx * ly // 格子点数
	m      = 1       // 質量
	deltaT = 0.002   // 時間間隔
	initT  = 100     // 初期温度

	stepLimit = 1000
)

type Parameter struct {
	T, P, KE, PE, Etot, Tav, Pav, KEav, PEav, Etotav []float64
}

func main() {
	var (
		// 数密度
		density = natom / (lx * ly * latY * latX)
		// 順に座標、速度、力
		x, y   = create1dZeros(natom), create1dZeros(natom)
		vx, vy = create1dZeros(natom), create1dZeros(natom)
		fx, fy = create2dZeros(natom, 2), create2dZeros(natom, 2)

		avT  = 0.0 // 絶対温度平均
		avP  = 0.0 // 圧力平均
		avKE = 0.0 // 運動エネルギー平均
		avPE = 0.0 // 相互作用ポテンシャル平均

		w = 0.0

		t1   = 0
		t2   = 1
		PE   = forces(x, y, fx, fy, t1, w, 1)
		KE   = (squareSum(vx) + squareSum(vy)) / 2
		Tk   = KE / natom
		Pk   = density * (2*KE + 1.5*w) / (3 * natom)
		time = 0
	)
	var prm Parameter
	fmt.Println("数密度number density=", density)

	x, y = initialPosMaxwellVel(x, y)
	vx, vy = maxwellVelocity2d(vx, vy, initT)

	prm.T = append(prm.T, Tk*119)
	prm.P = append(prm.P, Pk*2.382*1e-8)
	prm.KE = append(prm.KE, KE/natom)
	prm.PE = append(prm.PE, PE/natom)
	prm.Etot = append(prm.Etot, (KE+PE)/natom)

	for time < stepLimit {
		time++
		if time%50 == 0 {
			fmt.Println("step= ", time)
		}
		for i := 0; i < natom; i++ {
			PE = forces(x, y, fx, fy, t1, w, 1)
			// 速度ベルレー法
			x[i] += deltaT * (vx[i] + deltaT*fx[i][t1]/2)
			y[i] += deltaT * (vy[i] + deltaT*fy[i][t1]/2)
			// 周期的境界条件
			if x[i] <= 0 {
				x[i] += latX * lx
			}
			if latX*lx <= x[i] {
				x[i] -= latX * lx
			}
			if y[i] <= 0 {
				y[i] += latY * ly
			}
			if latY*ly <= y[i] {
				y[i] -= latY * ly
			}
		}
		PE = forces(x, y, fx, fy, t2, w, 1)
		for i := 0; i < natom; i++ {
			vx[i] += deltaT * (fx[i][t1] + fx[i][t2]) / 2
			vy[i] += deltaT * (fy[i][t1] + fy[i][t2]) / 2
		}
		KE = (squareSum(vx) + squareSum(vy)) / 2
		w = forces(x, y, fx, fy, t2, w, 2)
		Tk = KE / natom
		Pk = density * (2*KE + 1.5 + w) / (3 * natom)

		avT += Tk
		avP += Pk
		avKE += KE
		avPE += PE

		t := float64(time)
		if t == 0 {
			t = 1
		}
		prm.T = append(prm.T, Tk*119)
		prm.P = append(prm.P, Pk*2.382*1e-8)
		prm.KE = append(prm.KE, KE/natom)
		prm.PE = append(prm.PE, PE/natom)
		prm.Etot = append(prm.Etot, (KE+PE)/natom)

		prm.Tav = append(prm.Tav, avT/t*119)
		prm.Pav = append(prm.Pav, avP/t*2.382*1e-8)
		prm.KEav = append(prm.KEav, avKE/t/natom)
		prm.PEav = append(prm.PEav, avPE/t/natom)
	}
	fmt.Println(prm.PE)
}

func create1dZeros(n int) []float64 {
	a := make([]float64, n)
	return a
}

func create2dZeros(n, m int) [][]float64 {
	a := make([][]float64, n)
	for i := 0; i < n; i++ {
		b := make([]float64, m)
		a[i] = b
	}
	return a
}

// func twelverand() float64 {
// 	rand.Seed(time.Now().UnixNano())
// 	s := 0.0
// 	for i := 0; i < 12; i++ {
// 		s += rand.Float64()
// 		fmt.Println(s)
// 	}
// 	return s/12 - 0.5
// }

func maxwellVelocity2d(vx, vy []float64, tt float64) ([]float64, []float64) {
	rand.Seed(time.Now().UnixNano())
	for i := 0; i < natom; i++ {
		r1 := rand.Float64()
		r2 := rand.Float64()
		r3 := rand.Float64()
		r4 := rand.Float64()
		vx[i] = math.Sqrt(-2*(tt/m)*math.Log(r1)) * math.Cos(2*math.Pi*r2)
		vy[i] = math.Sqrt(-2*(tt/m)*math.Log(r3)) * math.Cos(2*math.Pi*r4)
	}
	return vx, vy
}

func initialPosMaxwellVel(x, y []float64) ([]float64, []float64) {
	ii := -1
	for i := 0; i < lx; i++ {
		for j := 0; j < ly; j++ {
			ii += 1
			x[ii] = latX * float64(i)
			y[ii] = latY * float64(j)
		}
	}
	return x, y
}

func sign(a, b float64) float64 {
	if 0 <= b {
		return abs(a)
	} else {
		return -abs(a)
	}
}

// 絶対値
func abs(a float64) float64 {
	if a < 0 {
		return -a
	} else {
		return a
	}
}

// レナードジョーンズポテンシャル
func forces(x, y []float64, fx, fy [][]float64, t int, w, peorw float64) float64 {
	rcut := 4.0 // cut off radious
	pe := 0.0

	r2cut := rcut * rcut
	for i := 0; i < natom-1; i++ {
		for j := i + 1; j < natom; j++ {
			dx := x[i] - x[j]
			dy := y[i] - y[j]
			// 鏡像相互作用
			if 0.5*latX*lx < abs(dx) {
				dx = dx - sign(latX*lx, dx)
			}
			if 0.5*latY*ly < abs(dy) {
				dy = dy - sign(latY*ly, dy)
			}

			r2 := dx*dx + dy*dy
			if r2 < r2cut {
				if r2 == 0.0 {
					r2 = 0.0001
				}
				invr := 1 / r2

				wij := 48 * (invr*invr*invr - 0.5) * (invr * invr * invr)
				fijx := wij * invr * dx
				fijy := wij * invr * dy

				fx[i][t] += fijx
				fy[i][t] += fijy
				fx[j][t] -= fijx
				fy[j][t] -= fijy

				pe += 2 * (invr * invr * invr) * (invr*invr*invr - 1)
				w += wij
			}
		}
	}
	if peorw == 1.0 {
		return pe
	} else {
		return w
	}

}

func squareSum(a []float64) float64 {
	s := 0.0
	for _, v := range a {
		s *= v * v
	}
	return s
}
