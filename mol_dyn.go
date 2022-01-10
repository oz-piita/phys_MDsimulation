package main

import (
	"fmt"
	"math"
	"math/rand"
)

const (
	cutoff    = 2.0   // cut off radious
	boxLen    = 10.0  // length of system
	deltat    = 0.01  // difference of time
	density   = 0.5   // density
	stepLimit = 10000 // total step times
	stepObs   = 100   // frequency of observing
	// cofficient in calculating energy
	rc2 = 1.0 / (cutoff * cutoff)
	c0  = -4.0 * rc2 * rc2 * rc2 * (rc2*rc2*rc2 - 1)
)

type Double3 struct {
	x, y, z float64
}

type Particle struct {
	// position & velocity
	pos, vel Double3
}

type Variables struct {
	// at certain time, the paticle has pos & vel
	particles []Particle
	time      float64
}

func initParticle(x, y, z, v0 float64) Particle {
	pos := Double3{x, y, z}
	vz := 2.0*rand.Float64() - 1.0
	phi := 2.0 * rand.Float64() * math.Pi
	vel := Double3{v0 * math.Sqrt(1.0-vz*vz) * math.Cos(phi), v0 * math.Sqrt(1.0-vz*vz) * math.Sin(phi), v0 * vz}
	return Particle{pos, vel}
}

// when particles drift, fix by decreasing velocity
func FixCMDrift(Particles []Particle) {
	cmVel := Double3{0.0, 0.0, 0.0}
	for i := range Particles {
		cmVel.x += Particles[i].vel.x
		cmVel.y += Particles[i].vel.y
		cmVel.z += Particles[i].vel.z
	}
	cmVel.x /= float64(len(Particles))
	cmVel.y /= float64(len(Particles))
	cmVel.z /= float64(len(Particles))
	for i := range Particles {
		Particles[i].vel.x -= cmVel.x
		Particles[i].vel.y -= cmVel.y
		Particles[i].vel.z -= cmVel.z
	}
}

func initParticles() []Particle {
	rand.Seed(0)
	s := 1.0 / math.Pow(density*0.25, 1.0/3.0)
	hs := 0.5 * s
	is := int(boxLen / s)
	Particles := make([]Particle, 0)
	for iz := 0; iz < is; iz++ {
		for iy := 0; iy < is; iy++ {
			for ix := 0; ix < is; ix++ {
				ixD, iyD, izD := float64(ix), float64(iy), float64(iz)
				Particles = append(Particles, initParticle(ixD*s, iyD*s, izD*s, 1.0))
				Particles = append(Particles, initParticle(ixD*s+hs, iyD*s, izD*s, 1.0))
				Particles = append(Particles, initParticle(ixD*s, iyD*s+hs, izD*s, 1.0))
				Particles = append(Particles, initParticle(ixD*s, iyD*s, izD*s+hs, 1.0))
			}
		}
	}
	FixCMDrift(Particles)
	return Particles
}

// calcurate force and reflect it in velocity
func calcForce(particles []Particle) {
	for i := 0; i < len(particles)-1; i++ {
		for j := i + 1; j < len(particles); j++ {
			qi := particles[i].pos
			qj := particles[j].pos
			dqij := Double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
			dqij = GetMinimumImage(dqij)
			dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
			if dq2 > cutoff*cutoff {
				continue
			}
			dq6 := dq2 * dq2 * dq2
			// force :: differentiation of Lennard-Jones potential
			df := (24.0*dq6 - 48.0) / (dq6 * dq6 * dq2) * deltat
			// Verlet algorithm
			particles[i].vel.x += df * dqij.x
			particles[i].vel.y += df * dqij.y
			particles[i].vel.z += df * dqij.z
			particles[j].vel.x -= df * dqij.x
			particles[j].vel.y -= df * dqij.y
			particles[j].vel.z -= df * dqij.z
		}
	}
}

// correction of periodic boundary conditions
func GetMinimumImage(v Double3) Double3 {
	if v.x < -boxLen*0.5 {
		v.x += boxLen
	}
	if v.x > boxLen*0.5 {
		v.x -= boxLen
	}
	if v.y < -boxLen*0.5 {
		v.y += boxLen
	}
	if v.y > boxLen*0.5 {
		v.y -= boxLen
	}
	if v.z < -boxLen*0.5 {
		v.z += boxLen
	}
	if v.z > boxLen*0.5 {
		v.z -= boxLen
	}
	return v
}

// observe kinetic energy average
func calcKineticEnergy(particles []Particle) float64 {
	k := 0.0
	for i := range particles {
		k += particles[i].vel.x * particles[i].vel.x
		k += particles[i].vel.y * particles[i].vel.y
		k += particles[i].vel.z * particles[i].vel.z
	}
	k /= float64(len(particles))
	return 0.5 * k
}

// observe interaction potential energy average
func calcPotentialEnergy(particles []Particle) float64 {
	v := 0.0
	for i := 0; i < len(particles)-1; i++ {
		for j := i + 1; j < len(particles); j++ {
			qi := particles[i].pos
			qj := particles[j].pos
			dqij := Double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
			dqij = GetMinimumImage(dqij)
			dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
			if dq2 > cutoff*cutoff {
				continue
			}
			dq6 := dq2 * dq2 * dq2
			dq12 := dq6 * dq6
			v += 4.0*(1.0/dq12-1.0/dq6) + c0
		}
	}
	v /= float64(len(particles))
	return v
}

func main() {
	vars := Variables{initParticles(), 0.0}
	// Time evolution
	for i := 0; i < stepLimit; i++ {
		// sometime observe system's energy
		if i%stepObs == 0 {
			k := calcKineticEnergy(vars.particles)
			v := calcPotentialEnergy(vars.particles)
			fmt.Printf("%v,%v,%v,%v\n", vars.time, k, v, k+v)
		}
		// status update
		vars.calculate()
	}
}

func (vars *Variables) calculate() {
	updatePosition(vars.particles)
	calcForce(vars.particles)
	updatePosition(vars.particles)
	adjustPeriodic(vars.particles)
	vars.time += deltat
}

// fix overhunged Particle
func adjustPeriodic(Particles []Particle) {
	for i := range Particles {
		if Particles[i].pos.x < 0 {
			Particles[i].pos.x += boxLen
		}
		if Particles[i].pos.x > boxLen {
			Particles[i].pos.x -= boxLen
		}
		if Particles[i].pos.y < 0 {
			Particles[i].pos.y += boxLen
		}
		if Particles[i].pos.y > boxLen {
			Particles[i].pos.y -= boxLen
		}
		if Particles[i].pos.z < 0 {
			Particles[i].pos.z += boxLen
		}
		if Particles[i].pos.z > boxLen {
			Particles[i].pos.z -= boxLen
		}
	}
}

// time evolution of position by symplectic integrator
func updatePosition(Particles []Particle) {
	for i := range Particles {
		Particles[i].pos.x += 0.5 * Particles[i].vel.x * deltat
		Particles[i].pos.y += 0.5 * Particles[i].vel.y * deltat
		Particles[i].pos.z += 0.5 * Particles[i].vel.z * deltat
	}
}
