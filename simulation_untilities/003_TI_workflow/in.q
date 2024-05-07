#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

log		log.q_${lambda}

##### simulation variables

variable        seed 	equal  4134
variable        T       equal  300.0		#K 
variable        p       equal  0.0 		#bar 
variable        damp_T  equal  100 		#K 
variable        damp_P  equal  1000.0		#bar 
variable        dT      equal  1.0 		#fs 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      lj/cut/coul/long 12.5
bond_style      harmonic
angle_style     hybrid cosine/periodic fourier
dihedral_style  harmonic
improper_style  fourier

dielectric      1.0
pair_modify     tail yes mix arithmetic
special_bonds   lj/coul 0.0 0.0 1.0
box tilt        large

read_data       ../data.${topo}_${cation}
kspace_style	ewald 1.0e-4

include 	../../MOF_q.parm
include 	../../${cation}_q.parm

group		fram type < 5
group		cation type > 4
group		center type 6

fix		SPRING fram spring/self 10 	# harmonic restraint to framework
fix		SPRING2 center spring/self 1

group           Fe type 2
variable	n_Fe equal count(Fe)
variable	fu equal v_n_Fe/4
print		"${fu}" file fu.dat

#------------------------------------------------------------------------------#
# equilibration
#------------------------------------------------------------------------------#

reset_timestep	0

variable 	delta equal ${lambda}
fix 		ADAPT all adapt 0 atom charge v_delta scale yes

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair
thermo 		5000
velocity 	all create ${T} ${seed}

fix		f1 all nve
fix 		f2 all langevin ${T} ${T} ${damp_T} ${seed} zero yes
compute 	c1 all temp/com
fix_modify	f2 temp c1

run		200000
unfix		ADAPT

#------------------------------------------------------------------------------#
# sampling
#------------------------------------------------------------------------------#

reset_timestep	0

variable	elec equal 2*(elong+ecoul)
variable	dudlam equal v_elec/${lambda}
fix		ELEC all ave/time 250 1 250 v_dudlam file q_${lambda}.lmp

## since compute group/group is expensive on the fly, require rerun to decompose the energy terms
## need to resolve the rerun issue

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair v_elec
thermo		250

dump            run_custom_data all custom 250 q_dump_all_${lambda} id mol type xu yu zu
dump_modify     run_custom_data sort id

run 		500000
write_data	data.restart_${topo}_${cation}_${lambda} nocoeff




