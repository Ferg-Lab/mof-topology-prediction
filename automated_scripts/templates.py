
em_template = '''
log		log.deposit

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}
kspace_style 	{kstyle}

dielectric      1.0
pair_modify     tail yes mix arithmetic
special_bonds   lj/coul 0.0 0.0 1.0
read_data       data.{mof}_new &
				extra/bond/per/atom 40 &
				extra/angle/per/atom 40 &
				extra/dihedral/per/atom 40 &
				extra/improper/per/atom 40 &
				extra/special/per/atom 40

include			emin_{mof}.param

'''


deposit_template = '''
#------------------------------------------------------------------------------#
# deposit cations and minimization
#------------------------------------------------------------------------------#
variable		a equal lx-2
variable 		b equal ly-2
variable 		c equal lz-2
variable		aa equal xy
variable 		bb equal xz
variable 		cc equal yz

region 			box prism 2 ${{a}} 2 ${{b}} 2 ${{c}} ${{aa}} ${{bb}} ${{cc}} 

variable		n_cation equal {n_cation}

reset_timestep  0
fix             f0 all box/relax aniso ${{p}} vmax 0.001
thermo          10
thermo_style    custom step lx ly lz press pxx pyy pzz
min_style       cg
minimize        1e-10 1e-10 100000 100000
unfix           f0

reset_timestep  0

##### create and deposit the cation

molecule 		{cation} {cation}.txt offset {Toff} {Boff} {Aoff} {Doff} {Ioff}
fix 			DEPOSIT all deposit {n_cation} 0 1 ${{seed}} mol {cation} region box near 2

run 			{step}

unfix 			DEPOSIT

'''


nvt_template = '''
##### NVT Dynamics

reset_timestep  0

velocity 		all create ${{T}} ${{T}}

fix		 		f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

thermo_style 	custom step pe ke etotal evdwl ecoul epair emol econserve density vol
thermo 	 	 	5000

run 		 	{equi_time}

unfix 		 	f1
unfix 		 	f2
'''

npt_template = '''
##### NPT Dynamics

reset_timestep  0

velocity 		all create ${{T}} ${{T}}

fix		 	f1 all npt temp ${{T}} ${{T}} ${{damp_T}} tri ${{p}} ${{p}} ${{damp_P}}
compute 		c1 all temp/com
fix_modify		f1 temp c1

thermo_style 	custom step pe ke etotal evdwl ecoul epair emol econserve density vol
thermo 	 	 	5000

run 		 	{equi_time}
unfix 		 	f1

'''

nvt_npt_template = '''
##### NVT/NPT Dynamics

reset_timestep  0

velocity 		all create ${{T}} ${{T}}

fix 		 	f1 all nve
fix 		 	f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify 		f2 temp c1

thermo_style 	custom step pe ke etotal evdwl ecoul epair emol econserve density vol
thermo 	 	 	5000

run 		 	{equi_time}

unfix 		 	f1
unfix 		 	f2

fix 		 	f1 all npt temp ${{T}} ${{T}} ${{damp_T}} tri ${{p}} ${{p}} ${{damp_P}}
fix_modify 		f1 temp c1
run 		 	{equi_time}

unfix 		 	f1

'''

output_template = '''
#------------------------------------------------------------------------------#
# final equilibration and save data file
#------------------------------------------------------------------------------#
### save volume
reset_timestep  0
min_style 		cg
minimize		1e-15 1e-15 100000 100000

reset_timestep  0
fix             f0 all box/relax tri ${{p}} vmax 0.001
thermo          10
thermo_style    custom step lx ly lz press pxx pyy pzz
min_style       cg
minimize        1e-15 1e-15 100000 100000
unfix           f0

min_style 		cg
minimize		1e-15 1e-15 100000 100000

variable 		volume equal vol
print 		 	"The volume of {mof}: ${{volume}}" file volume_{mof}.dat

write_data 		data.emin_{mof}_{cation} pair ij

### For strain energy
reset_timestep  0


fix				f1 all nve
fix 			f2 all langevin 0 0 ${{damp_T}} ${{seed}}

thermo_style 	custom step pe ke etotal evdwl ecoul epair emol econserve density vol
thermo 	 	 	5000

variable		strain equal pe
fix 			STRAIN all ave/time 1 100 100 v_strain file STRAIN_{mof}.dat

run 			100

unfix			f1
unfix			f2
unfix			STRAIN

### For potential energy
reset_timestep	0
velocity 		all create ${{T}} ${{T}}
fix 		 	f1 all nve
fix 		 	f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
fix_modify 		f2 temp c1
thermo_style 	custom step pe ke etotal evdwl ecoul epair emol econserve density vol
thermo 	 	 	5000
variable 	 	potential equal pe
fix				PE all ave/time 1 50000 50000 v_potential file PE_{mof}.dat

run				50000

unfix			f1
unfix			f2
unfix			PE


'''


bonded_part1_template = '''
log		./{mof}/bonded/log.bonded_${{lambda}}

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}

pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	data.emin_{mof}_{job}
include 	 	bonded_{mof}.param

group			fram type < {fram_atoms}

fix				SPRING fram spring/self 10 	# harmonic restraint to framework

'''

bonded_part2_template = '''

#------------------------------------------------------------------------------#
# equilibration
#------------------------------------------------------------------------------#

reset_timestep	0

thermo_style	custom step etotal ke pe temp press density vol ebond eangle edihed eimp emol 
thermo 			5000

velocity 		all create ${{T}} ${{seed}}

fix				f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

run				200000

#------------------------------------------------------------------------------#
# sampling
#------------------------------------------------------------------------------#

reset_timestep	0

compute 		MOF_pe_a fram pe/atom bond angle dihedral improper	# per atom quantity
compute			MOF_pe all reduce sum c_MOF_pe_a	# reduce to bulk quantity

variable 		dudlam	equal c_MOF_pe/${{lambda}}

fix				BOND all ave/time 250 1 250 v_dudlam file ./{mof}/bonded/bond_${{lambda}}.lmp

thermo_style	custom step etotal ke pe temp press density vol c_MOF_pe v_dudlam 
thermo 			250

run 			500000

unfix 			f1
unfix 			f2

'''


lj_part1_template = '''
log		./{mof}/lj/log.lj_${{lambda}}

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}

pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	data.emin_{mof}_{job}
include 	 	non_bonded_{mof}_lj.param

group			fram type < {fram_atoms}

fix				SPRING fram spring/self 10 	# harmonic restraint to framework

'''

lj_part2_template = '''
#------------------------------------------------------------------------------#
# equilibration
#------------------------------------------------------------------------------#

reset_timestep	0

variable 		delta equal ${{lambda}}
fix 			ADAPT all adapt 0 pair lj/cut/soft lambda * * v_delta

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair
thermo 			5000
velocity 		all create ${{T}} ${{seed}}

fix				f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

run				500000
unfix			ADAPT

#------------------------------------------------------------------------------#
# sampling
#------------------------------------------------------------------------------#

reset_timestep	0

variable		ddlambda equal 0.001
variable		ddlambda_n equal -0.001

compute 		ALL1 all fep ${{T}} pair lj/cut/soft lambda * * v_ddlambda
compute 		ALL2 all fep ${{T}} pair lj/cut/soft lambda * * v_ddlambda_n
variable 		dudlam equal (c_ALL1[1]-c_ALL2[1])/0.002

fix				ALL all ave/time 250 1 250 v_dudlam file ./{mof}/lj/LJ_${{lambda}}.lmp

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair v_dudlam
thermo			250

run 			800000

unfix			f1
unfix			f2

'''

q_part1_template = '''
log		./{mof}/q/log.q_${{lambda}}

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}
kspace_style 	{kstyle}

dielectric 	 	1.0
pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	data.emin_{mof}_{job}
include 	 	non_bonded_{mof}.param

group			fram type < {fram_atoms}

fix				SPRING fram spring/self 10 	# harmonic restraint to framework

'''


q_part2_template = '''
#------------------------------------------------------------------------------#
# equilibration
#------------------------------------------------------------------------------#

reset_timestep	0

variable 		delta equal ${{lambda}}
fix 			ADAPT all adapt 0 atom charge v_delta scale yes

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair
thermo 			5000
velocity 		all create ${{T}} ${{seed}}

fix				f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

run				200000
unfix			ADAPT

#------------------------------------------------------------------------------#
# sampling
#------------------------------------------------------------------------------#

reset_timestep	0

variable		elec equal 2*(elong+ecoul)
variable		dudlam equal v_elec/${{lambda}}
fix				ELEC all ave/time 250 1 250 v_dudlam file ./{mof}/q/q_${{lambda}}.lmp

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair v_elec
thermo			250

run 			500000

unfix			f1
unfix			f2

'''


hr_with_mol_template = '''
log		./{mof}/hr/log.hr

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}
kspace_style 	{kstyle}

dielectric 	 	1.0
pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	data.emin_{mof}_{job}
include 	 	non_bonded_{mof}.param

group			fram type < {fram_atoms}
group			cation type > {cation_atoms}
group			center type {center_atom}

#------------------------------------------------------------------------------#
# ti/spring definition
#------------------------------------------------------------------------------#

reset_timestep	0

variable		t_s equal 1000000
variable		t_eq equal 500000

variable        t_eq equal ${{t_eq}}-1
variable        t_s equal ${{t_s}}+1

fix             fs1 fram ti/spring 10 ${{t_s}} ${{t_eq}} function 2
fix             fs2 center ti/spring 1 ${{t_s}} ${{t_eq}} function 2

variable        dE equal f_fs1+f_fs2
variable        lambda equal f_fs1[1]

fix				f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair f_fs1 f_fs2
thermo 			1000

#------------------------------------------------------------------------------#
# switching
#------------------------------------------------------------------------------#

### forward integration, save to forward.dat

run             ${{t_eq}}
fix             f4 all print 1 "${{dE}} ${{lambda}}" title "# dE lambda" &
                screen no file ./{mof}/hr/forward.dat
run             ${{t_s}}
unfix           f4

### backward integration, save to backward.dat

run             ${{t_eq}}
fix             f4 all print 1 "${{dE}} ${{lambda}}" title "# dE lambda" &
                screen no file ./{mof}/hr/backward.dat
run             ${{t_s}}
unfix           f4

print			"done" file ./{mof}/hr/done.dat

unfix			f1
unfix			f2

'''


hr_no_mol_template = '''
log		./{mof}/hr/log.hr

#------------------------------------------------------------------------------#
# simulation settings
#------------------------------------------------------------------------------#

##### simulation variables

variable        seed 	equal  {seed}
variable        T       equal  {temp}	    #K 
variable        p       equal  {pressure} 	#bar 
variable        damp_T  equal  100 		  	#K 
variable        damp_P  equal  1000.0	  	#bar 
variable        dT      equal  1.0 		 	#fs 
variable        kb      equal  0.0019872041 

## system setup ## 

units           real
atom_style      full
boundary        p p p

pair_style      {pstyle}
bond_style      {bstyle}
angle_style     {astyle}
dihedral_style  {dstyle}
improper_style  {istyle}
kspace_style 	{kstyle}

dielectric 	 	1.0
pair_modify 	tail yes mix arithmetic
special_bonds 	lj/coul 0.0 0.0 1.0

read_data 	 	data.emin_{mof}_{job}
include 	 	non_bonded_{mof}.param

group			fram type < {fram_atoms}

#------------------------------------------------------------------------------#
# ti/spring definition
#------------------------------------------------------------------------------#

reset_timestep	0

variable		t_s equal 1000000
variable		t_eq equal 500000

variable        t_eq equal ${{t_eq}}-1
variable        t_s equal ${{t_s}}+1

fix             fs1 fram ti/spring 10 ${{t_s}} ${{t_eq}} function 2

variable        dE equal f_fs1
variable        lambda equal f_fs1[1]

fix				f1 all nve
fix 			f2 all langevin ${{T}} ${{T}} ${{damp_T}} ${{seed}} zero yes
compute 		c1 all temp/com
fix_modify		f2 temp c1

thermo_style	custom step etotal ke pe temp press density vol evdwl etail ecoul elong epair f_fs1
thermo 			1000

#------------------------------------------------------------------------------#
# switching
#------------------------------------------------------------------------------#

### forward integration, save to forward.dat

run             ${{t_eq}}
fix             f4 all print 1 "${{dE}} ${{lambda}}" title "# dE lambda" &
                screen no file ./{mof}/hr/forward.dat
run             ${{t_s}}
unfix           f4

### backward integration, save to backward.dat

run             ${{t_eq}}
fix             f4 all print 1 "${{dE}} ${{lambda}}" title "# dE lambda" &
                screen no file ./{mof}/hr/backward.dat
run             ${{t_s}}
unfix           f4

print			"done" file ./{mof}/hr/done.dat

unfix			f1
unfix			f2

'''


