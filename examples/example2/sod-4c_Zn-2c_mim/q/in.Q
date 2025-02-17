
log		log.q_${lambda}

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  20902
variable        T       equal  300	    #K 
variable        p       equal  0.0 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      lj/cut/coul/long 12.500
bond_style      harmonic
angle_style     hybrid cosine/periodic fourier
dihedral_style  harmonic
improper_style  fourier
kspace_style 	ewald 0.000001

dielectric 	 	1.0
pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	../data.emin_sod-4c_Zn-2c_mim_fe
include 	 	../non_bonded_sod-4c_Zn-2c_mim.param

group			fram type < 6

fix				SPRING fram spring/self 10 	# harmonic restraint to framework


#------------------------------------------------------------------------------#
# equilibration
#------------------------------------------------------------------------------#

reset_timestep	0

variable 		delta equal ${lambda}
fix 			ADAPT all adapt 0 atom charge v_delta scale yes

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair
thermo 			5000
velocity 		all create ${T} ${seed}

fix				f1 all nve
fix 			f2 all langevin ${T} ${T} ${damp_T} ${seed} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

run				200000
unfix			ADAPT

#------------------------------------------------------------------------------#
# sampling
#------------------------------------------------------------------------------#

reset_timestep	0

variable		elec equal 2*(elong+ecoul)
variable		dudlam equal v_elec/${lambda}
fix				ELEC all ave/time 250 1 250 v_dudlam file q_${lambda}.lmp

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair v_elec
thermo			250

run 			500000

unfix			f1
unfix			f2

