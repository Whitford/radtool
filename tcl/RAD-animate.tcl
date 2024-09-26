##########################################################################
# The RADtool software is a plugin for VMD, which is designed for 
# the analysis of ribosome structures.
#   Copyright (C) 2021  Paul C. Whitford
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
##########################################################################

# For questions, contact Paul Whitford (p.whitford@northeastern.edu)

namespace eval RTANIMATE {
 variable pi
 set pi 3.14159265358979323846
 variable rtd
 variable seqrtt
 variable deg2rad
 variable rad2deg
 variable headst
 set deg2rad [expr $pi/180.0]
 set rad2deg [expr 1.0/$deg2rad]

 proc RT-animate {optin ID_SSUrRNA ID_REF_SrRNA_1 ID_REF_SrRNA_2 vizseq Nrot_head Ntilt_head Ntrans_head Nrot_body Ntilt_body Ntrans_body } {
  variable seqrtt
  set seqrtt $vizseq 
  upvar $optin options
  # load in atom definitions for calculating angles
  foreach i { body0 body1 body2 head0 head1 head2 body_tilt0_0 body_tilt0_1 head_tilt0_0 head_tilt0_1 body_d head_d } { 
   set $i $::RADTOOL::RTD($i)
  }

 
  # the ID refers to the reference SSU rRNA structure.  We will separately move the full subunit (for body rotation)
  # and just the head
  set body [atomselect $ID_SSUrRNA "all"]
  set head [atomselect $ID_SSUrRNA "resid 928 to 1389"]

 # grab info about last model 
  set RTHIST $::RADTOOL::radenv(RADtoolhistory)
  set endmol [lindex $RTHIST end]

  lassign $endmol rot_body_end tilt_body_end tiltdir_body_end trans_body_end bodyorigin rot_head_end tilt_head_end tiltdir_head_end trans_head_end headorigin hf hl modelID_body_end modelID_head_end smallchain_end LSUID_end SSUID_end 

  if { ! [info exists options(ONLYSSU)] } {
   ############# DO FULL SSU FIRST ########################################
   puts "\n\nAnimation Information: Body"
   set zero_dir [ ::RADTOOL::get_one_vector $ID_SSUrRNA $body_tilt0_0 $body_tilt0_1 ]
   animatestage $ID_SSUrRNA $body $body0 $body1 $body2 $rot_body_end $tilt_body_end $tiltdir_body_end $trans_body_end $bodyorigin $Nrot_body $Ntilt_body $Ntrans_body $zero_dir 
  }
 ############# DO HEAD ########################################
  puts "\n\nAnimation Information: Head"
  set zero_dir [ ::RADTOOL::get_one_vector $ID_SSUrRNA $head_tilt0_0 $head_tilt0_1 ]
  animatestage $ID_SSUrRNA $head $head0 $head1 $head2 $rot_head_end $tilt_head_end $tiltdir_head_end $trans_head_end $headorigin $Nrot_head $Ntilt_head $Ntrans_head $zero_dir
  animate delete beg 0 end 0 $ID_SSUrRNA

  checkfinalorientation [info exists options(ONLYSSU)] $ID_REF_SrRNA_1 $ID_REF_SrRNA_2 $ID_SSUrRNA 928 1389
  animate speed 0.9
  animate goto 0
  animate style Rock
  animate forward
  mol on $ID_SSUrRNA
  mol top $ID_SSUrRNA

 }
 
 proc RT-animate2 {args} {
  variable rtd
  variable seqrtt
  array unset rtd  
  if { [llength $args] == 5 } {
   lassign $args startmol endmol Nsteps modn seqrtt
  } elseif { [llength $args] == 0} {
   #only allow default animation if called from the command line
   puts stderr "Going to animate between two models"
   puts stderr "What is the ID of the starting model? Leave blank to animate two most recent models."
   set startmol [gets stdin]
   if { $startmol eq "" } {
    set startmol -1
   }
   if { $startmol != -1 } {
    puts stderr "What is the ID of the ending model?"
    set endmol [gets stdin]
   } else {
    set endmol -1
   } 

   puts stderr "How many frames for each step? \[10\]"
   set Nsteps [gets stdin]
   if { $Nsteps eq "" } {
    set Nsteps 10
   }
   puts stderr "Which structure should be used for animation?\n0=classical\n1=starting model"
   set modn [gets stdin]
   if { $modn eq "" } {
    set modn 0
   }
   puts stderr "Which sequence should be used for animation?\n\t0=translate, rotate, tilt\n\t1=translate, tilt, rotate\n\t2=rotate, tilt, translate\n\t3=tilt, rotate, translate"
   set seqs [gets stdin]
   if { $seqs eq "" } {
    set seqrtt "r_t_tr"
   } elseif {$seqs eq "0"} {
    set seqrtt "tr_r_t"
   } elseif {$seqs eq "1"} {
    set seqrtt "tr_t_r"
   } elseif {$seqs eq "2"} {
    set seqrtt "r_t_tr"
   } elseif {$seqs eq "3"} {
    set seqrtt "t_r_tr"
   } else {
    error "Invalid value: \"$seqs\"\n\tSupported values: r_t_rt, t_r_tr, tr_r_t, tr_t_r"
   }

  } else {
   error "RT-animate2 can only be called with 0, or 5, values."
  } 

  if { ![string is integer "$startmol"] } {
   error "Model ID of starting model must be an integer (mol #). Found $startmol"
  }

  if { ![string is integer "$endmol"] } {
   error "Model ID of ending model must be an integer (mol #). Found $endmol"
  }
 
  if { ![string is integer "$Nsteps"] } {
   error "Number of steps must be an integer. Found $Nsteps"
  }
  if { $Nsteps < 1 } {
   error "Number of steps must be positive. Found $Nsteps"
  }

  set RTHIST $::RADTOOL::radenv(RADtoolhistory)
  set RTHISL [llength $RTHIST]
  if { $startmol == -1 && $endmol ==-1} {
   # by default, just animate from the last 2
   if { $RTHISL < 2 } {
    error "When using -animate2, there must be at least 2 models in the history"
   }
   set endmol [lindex $RTHIST end]
   set startmol [lindex $RTHIST end-1]
  } elseif { $startmol == -1 || $endmol == -1 } {
   error "When using -animate2, must give two model indices, or none"
  } else {
   # must both be set
   if { $endmol >= $RTHISL } {
    error "-animate2 endpoint index $endmol does not exist"
   }
   if { $startmol >= $RTHISL } {
    error "-animate2 startpoint index $startmol does not exist"
   }
   set endmol [lindex $RTHIST $endmol]
   set startmol [lindex $RTHIST $startmol]
  }
  # retrieve info about earlier models
  set startmol [join $startmol]
  lassign $startmol rot_body_begin tilt_body_begin tiltdir_body_begin trans_body_begin bodyorigin_begin rot_head_begin tilt_head_begin tiltdir_head_begin trans_head_begin headorigin_begin hf hl modelID_body_begin modelID_head_begin smallchain_begin LSUID_begin SSUID_begin 
 
  if {$LSUID_begin == -100} {
   # must be a head-only case
   set ONLYSSU 0
  } else {
   set ONLYSSU 1
  }
  set endmol [join $endmol]
  lassign $endmol rot_body_end tilt_body_end tiltdir_body_end trans_body_end bodyorigin_end rot_head_end tilt_head_end tiltdir_head_end trans_head_end headorigin_end hft hlt modelID_body_end modelID_head_end smallchain_end LSUID_end SSUID_end
  if { ! $ONLYSSU } { 
   set ERbody [ERfromEuler $rot_body_begin $tilt_body_begin $tiltdir_body_begin $rot_body_end $tilt_body_end $tiltdir_body_end] 
   puts "Euler-Rodrigues angle formed by the two body orientations is [format %.2f $ERbody] degrees"
  }
  set ERhead [ERfromEuler $rot_head_begin $tilt_head_begin $tiltdir_head_begin $rot_head_end $tilt_head_end $tiltdir_head_end] 
  puts "Euler-Rodrigues angle formed by the two head orientations (relative to the same body orientation) is [format %.2f $ERhead] degrees"
  # left this way so it may be customized later
  set Nrot_body   $Nsteps
  set Ntilt_body  $Nsteps
  set Ntiltd_body $Nsteps 
  set Ntrans_body $Nsteps
  set Nrot_head   $Nsteps
  set Ntilt_head  $Nsteps
  set Ntiltd_head $Nsteps  
  set Ntrans_head $Nsteps
  if {$ONLYSSU==0} {
   # since only SSU, add one to head rot display, so that initial frame is shown
   incr Nrot_head
  } else {
   incr Nrot_body
  }
  # load the reference model and then return the origins and centers back to the reference frame
  set ID_SSUrRNA  [mol new "$::RADTOOL::radenv(ROTATIONPATH)/../share/reference_models/4v9d.16S.pdb" waitfor -1] 
  set selbegin [atomselect $modelID_body_begin "not resid 928 to 1389"]
  set selend [atomselect $modelID_body_end "not resid 928 to 1389"]
  set selref [atomselect $ID_SSUrRNA "not resid 928 to 1389"]
  set M [measure fit $selbegin $selref]
  set headorigin_begin [coordtrans $M $headorigin_begin]
  set trans_head_begin [vectrans $M $trans_head_begin]
  set M2 [measure fit $selend $selref]
  set headorigin_end [coordtrans $M2 $headorigin_end]
  set trans_head_end [vectrans $M2 $trans_head_end]
  $selbegin delete
  $selend delete
  $selref delete
 
  # copy the initial SSU model and reorient to the same reference
  set ID_SSUrRNA2 [::RADTOOL::restorenumbering $SSUID_begin 0]
  set op(OUTPUT) "stderr"
  ::RADTOOL::findlargestblock op $ID_SSUrRNA2 $smallchain_begin 0
  set sel [atomselect $ID_SSUrRNA2 "nucleic and chain \"$smallchain_begin\" and not resid $hf to $hl"]
  $sel move $M
  $sel delete
  set selheadbegin [atomselect $modelID_head_begin "resid 928 to 1389"]
  set selheadref [atomselect $ID_SSUrRNA "resid 928 to 1389"]
  set M [measure fit $selheadbegin $selheadref]
  set sel [atomselect $ID_SSUrRNA2 "nucleic and chain \"$smallchain_begin\" and resid $hf to $hl"]
  $sel move $M
  $sel delete
  $selheadbegin delete
  $selheadref delete
  
  # load in atom definitions for calculating angles
  foreach i { body0 body1 body2 head0 head1 head2 body_tilt0_0 body_tilt0_1 head_tilt0_0 head_tilt0_1 body_d head_d } {
   set $i $::RADTOOL::RTD($i)
  }
 
  # the ID refers to the reference SSU rRNA structure.  We will separately move the full subunit (for body rotation)
  # and just the head
  set body [atomselect $ID_SSUrRNA "all"]
  set head [atomselect $ID_SSUrRNA "resid 928 to 1389"]
  set bodymod [atomselect $ID_SSUrRNA2 "nucleic and chain \"$smallchain_begin\" "]
  set headmod [atomselect $ID_SSUrRNA2 "nucleic and chain \"$smallchain_begin\" and resid $hf to $hl"]
  if {$ONLYSSU == 1} { 
   puts "Will try to display a rigid-body representation of the 
 structural differences in the last two models analyzed with RADtool. 
 To simplify the display, the SSU and LSU will be redrawn. The 
 starting conformation (Mol ID $LSUID_begin and $SSUID_begin) is cyan.
 The end conformation (Mol ID $LSUID_end and $SSUID_end) is yellow."
  } else {
   puts "Will try to display a rigid-body representation of the 
 structural differences in the last two models analyzed with RADtool. 
 To simplify the display, the SSU will be redrawn. The 
 starting conformation (Mol ID $SSUID_begin) is cyan.
 The end conformation (Mol ID $SSUID_end) is yellow."
  }
  animateall $ONLYSSU $modn $ID_SSUrRNA2 $bodymod $headmod $ID_SSUrRNA $body $head $body0 $body1 $body2 $rot_body_begin $tilt_body_begin $tiltdir_body_begin $trans_body_begin $rot_body_end $tilt_body_end $tiltdir_body_end $trans_body_end $bodyorigin_begin $bodyorigin_end $Nrot_body $Ntilt_body $Ntiltd_body $Ntrans_body $body_tilt0_0 $body_tilt0_1 $head0 $head1 $head2 $rot_head_begin $tilt_head_begin $tiltdir_head_begin $trans_head_begin $rot_head_end $tilt_head_end $tiltdir_head_end $trans_head_end $headorigin_begin $headorigin_end  $Nrot_head $Ntilt_head $Ntiltd_head $Ntrans_head $head_tilt0_0 $head_tilt0_1
  if {$modn == 0 } {
   checkfinalorientation $ONLYSSU $modelID_body_end $modelID_head_end $ID_SSUrRNA 928 1389
   # need to add a check for when we are animating the cores
  } 
  # just show tubes for the two structures' LSU and SSU rRNA
  mol off all
  if {$ONLYSSU == 1} { 
   mol on $LSUID_begin
   mol on $LSUID_end
   mol showrep $LSUID_begin 0 0
   mol showrep $LSUID_end   0 0
   mol showrep $LSUID_begin 1 1
   mol showrep $LSUID_end   1 1
   mol modcolor 1 $LSUID_begin ColorID 10
   mol modcolor 1 $LSUID_end ColorID 4
  }
 
  # show the endpoint configurations 
  mol on $SSUID_begin
  mol on $SSUID_end
  # hide lines
  mol showrep $SSUID_begin 0 0
  mol showrep $SSUID_end   0 0
  # show tubes
  mol showrep $SSUID_begin 1 1
  mol showrep $SSUID_end   1 1
  # assign SSU colors
  mol modcolor 1 $SSUID_begin ColorID 10
  mol modcolor 1 $SSUID_end ColorID 4 

  if {$modn ==0} { 
   # rigid-body animation
   animate delete beg 0 end 0 $ID_SSUrRNA
   mol rename $ID_SSUrRNA "animation"
   ::RADTOOL::drawtube $ID_SSUrRNA "928 to 1389" 1.4 15
   ::RADTOOL::drawreptube $ID_SSUrRNA "0 to 927 1390 to 10000" 1.4 1 "on"
   mol on $ID_SSUrRNA
  } else {
   # animation based on beginning model
   animate delete beg 0 end 0 $ID_SSUrRNA2
   mol rename $ID_SSUrRNA2 "animation"
   ::RADTOOL::drawtube $ID_SSUrRNA2 "nucleic and chain $smallchain_begin and resid $hf to $hl" 1.4 15
   ::RADTOOL::drawreptube $ID_SSUrRNA2 "nucleic and chain $smallchain_begin and not resid $hf to $hl" 1.4 1 "on"
   mol on $ID_SSUrRNA2
  }
 
  animate speed 0.9
  animate goto 0
  animate style Rock
  animate forward
  ::RADTOOL::orientribosome
 }
 
 proc animatestage {MOLID SELID res0 res1 res2 rot_end tilt_end tiltdir_end translate_end origin Nrot Ntilt Ntrans zero_dir } { 
 
  variable pi
  variable rtd
  variable seqrtt

  array unset rtd  
  set fnstart [expr [molinfo $MOLID get numframes] - 1 ]
  set fnn  $fnstart

  if {$Nrot < 0} {
   puts "Number of rotation values can not be negative."
   return
  } elseif {$Nrot > 0} {
   set stepr [expr 1.00*$rot_end/$Nrot]
  } else {
   set stepr 0
   set Nrot 0
  }
 
  if { ($tilt_end > 180) || ($tilt_end < 0)} {
   puts "Tilt angle is only defined for 0 to 180 degrees"
   return
  }
  if {$Ntilt < 0} {
   puts "Number of tilt values can not be negative."
   return
  } elseif {$Ntilt > 0} {
   set stept [expr 1.00*$tilt_end/$Ntilt]
  } else {
   set stept 0
   set Ntilt 0
  }

  set tiltdir_end [expr $pi*$tiltdir_end/180.0 ]
  set translate $translate_end 
  if {$Ntrans < 0} {
   error "Number of translation points can not be negative."
  } elseif {$Ntrans > 0} {
   set transstep [vecscale [expr 1.0/$Ntrans] $translate ]
  } else {
   set Ntrans 0
  }
  # make a sequence of values for later
  set rtd($fnn,rot) 0
  set rtd($fnn,tilt) 0
  set rtd($fnn,tiltd) $tiltdir_end
  set rtd($fnn,trans) {0 0 0}
  if { $seqrtt eq "r_t_tr"} {
   set fnn [genr $fnn $stepr $Nrot] 
   set fnn [gent $fnn $stept $Ntilt] 
   set fnn [gentr $fnn $transstep $Ntrans] 
  } elseif { $seqrtt eq "t_r_tr"} {
   set fnn [gent $fnn $stept $Ntilt] 
   set fnn [genr $fnn $stepr $Nrot] 
   set fnn [gentr $fnn $transstep $Ntrans] 
  } elseif { $seqrtt eq "tr_r_t"} {
   set fnn [gentr $fnn $transstep $Ntrans] 
   set fnn [genr $fnn $stepr $Nrot] 
   set fnn [gent $fnn $stept $Ntilt] 
  } elseif { $seqrtt eq "tr_t_r"} {
   set fnn [gentr $fnn $transstep $Ntrans] 
   set fnn [gent $fnn $stept $Ntilt] 
   set fnn [genr $fnn $stepr $Nrot] 
  } else {
   error "Invalid value given with -animate: \"$seqrtt\"\n\tSupported values: r_t_rt, t_r_tr, tr_r_t, tr_t_r"
  }
 
  # start by rotating and tilting 
  lassign [::RADTOOL::get_vectors $MOLID $res0 $res1 $res2] V3 V1 V2 
  
  # these are the rotations based on the e coli structure
 
  set zero_tilt_dir $zero_dir
  set newone [vecnorm [veccross $V1 $zero_tilt_dir]]
  set newzero [vecnorm  [veccross $newone $V1]]
  
  # rotate about the axis of rotation.  
  set rot_vec $V1
  set vec_end_rot [vecadd $origin $rot_vec]
  $SELID frame end
  for {set ir $fnstart} {$ir <= $fnn} {incr ir} {
   animate dup frame $fnstart $MOLID
   $SELID frame last
   moveelement $SELID $origin $vec_end_rot $rtd($ir,rot) $rtd($ir,tiltd) $rtd($ir,tilt) $rtd($ir,trans) $newone $newzero  
  }
 }

 proc genr {fnn stepr Nrot} {
  variable rtd
  set pframe $fnn
  set rb $rtd($pframe,rot)
  if { $Nrot > 0} {
   puts "    Rotation begins with frame [expr $fnn+1]"
   for {set ir 1} {$ir < $Nrot+1} {incr ir} { 
    incr fnn
    set rtd($fnn,rot)   [expr $rb+$stepr*$ir ]
    set rtd($fnn,tilt)  $rtd($pframe,tilt)
    set rtd($fnn,tiltd) $rtd($pframe,tiltd)
    set rtd($fnn,trans) $rtd($pframe,trans)
   }
  }
  return $fnn
 }

 proc gent {fnn stept Ntilt} {
  variable rtd
  set pframe $fnn
  set tb $rtd($pframe,tilt)

  if {$Ntilt > 0 } {
   puts "    Tilting begins with frame [expr $fnn+1]"

   for {set ir 1} {$ir < $Ntilt+1} {incr ir} { 
    incr fnn
    set rtd($fnn,rot) $rtd($pframe,rot) 
    set rtd($fnn,tilt) [expr $tb+$stept*$ir ]
    set rtd($fnn,tiltd) $rtd($pframe,tiltd) 
    set rtd($fnn,trans) $rtd($pframe,trans) 
   }
  }
  return $fnn 
 }

 proc gentr {fnn transstep Ntrans} {
  variable rtd
  set pframe $fnn
  set trb $rtd($pframe,trans)
  if {$Ntrans >0} {
   puts "    Translation starts at frame [expr $fnn+1]"
  }
  for {set ir 1} {$ir < $Ntrans+1} {incr ir} { 
   incr fnn
   set rtd($fnn,rot)   $rtd($pframe,rot)   
   set rtd($fnn,tilt)  $rtd($pframe,tilt)  
   set rtd($fnn,tiltd) $rtd($pframe,tiltd)  
   set rtd($fnn,trans) [vecadd $trb [vecscale $ir $transstep]]
  }
  return $fnn
 }
 
 proc animateall {ONLYSSU modn MOLIDMOD SELIDALLMOD SELIDHEADMOD MOLID SELIDALL SELIDHEAD b_res0 b_res1 b_res2 b_rot_begin b_tilt_begin b_tiltdir_begin b_translate_begin b_rot_end b_tilt_end b_tiltdir_end b_translate_end b_origin_begin b_origin_end b_Nrot b_Ntilt b_Ntiltd b_Ntrans body_tilt0_0 body_tilt0_1 h_res0 h_res1 h_res2 h_rot_begin h_tilt_begin h_tiltdir_begin h_translate_begin h_rot_end h_tilt_end h_tiltdir_end h_translate_end h_origin_begin h_origin_end h_Nrot h_Ntilt h_Ntiltd h_Ntrans head_tilt0_0 head_tilt0_1 } { 
  ############# SET HEAD, RELATIVE TO BODY ###############
  variable rtd 
  variable pi
  variable deg2rad
  variable rad2deg
  variable headst
  variable seqrtt
  # it is assumed that our starting frame is always the last one in the MOLID
  if {$b_tilt_begin == 0} {
   set Ntiltd 0
   # nothing to show during the tilt direction change.  So, don't increment.
  }
 
  set b_stepr [getrotstep $b_rot_begin $b_rot_end $b_Nrot]
  set b_stept [gettiltstep $b_tilt_begin $b_tilt_end $b_Ntilt]
  lassign [getoriginstep $b_tilt_begin $b_tilt_end $b_origin_begin $b_origin_end $b_Ntilt] b_origin_begin b_origin_end b_origin_step
  lassign [gettiltdstep $b_tiltdir_begin $b_tiltdir_end $b_Ntiltd] b_tiltdir_begin b_tiltdir_end b_steptd
  set b_transstep [gettransstep $b_translate_begin $b_translate_end $b_Ntrans]
 
  set h_stepr [getrotstep $h_rot_begin $h_rot_end $h_Nrot]
  set h_stept [gettiltstep $h_tilt_begin $h_tilt_end $h_Ntilt]
  lassign [getoriginstep $h_tilt_begin $h_tilt_end $h_origin_begin $h_origin_end $h_Ntilt] h_origin_begin h_origin_end h_origin_step
  lassign [gettiltdstep $h_tiltdir_begin $h_tiltdir_end $h_Ntiltd] h_tiltdir_begin h_tiltdir_end h_steptd 
  set h_transstep [gettransstep $h_translate_begin $h_translate_end $h_Ntrans]
 
  set headst 0 

  initrtd $ONLYSSU $b_rot_begin $b_tilt_begin $b_tiltdir_begin $b_translate_begin $b_origin_begin $h_rot_begin $h_tilt_begin $h_tiltdir_begin $h_translate_begin $h_origin_begin 

  if {$ONLYSSU == 1} { 
  # make a sequence of values for later 
   incr fnn
   if { $seqrtt eq "r_t_tr"} {
    set fnn [genbr $fnn $b_stepr $b_Nrot]
    set fnn [genbt $fnn $b_Ntilt $b_tilt_end $b_tiltdir_end $b_origin_step]
    set fnn [genbtrans $fnn $b_transstep $b_Ntrans]
   } elseif { $seqrtt eq "t_r_tr"} {
    set fnn [genbt $fnn $b_Ntilt $b_tilt_end $b_tiltdir_end $b_origin_step]
    set fnn [genbr $fnn $b_stepr $b_Nrot]
    set fnn [genbtrans $fnn $b_transstep $b_Ntrans]
   } elseif { $seqrtt eq "tr_r_t"} {
    set fnn [genbtrans $fnn $b_transstep $b_Ntrans]
    set fnn [genbr $fnn $b_stepr $b_Nrot]
    set fnn [genbt $fnn $b_Ntilt $b_tilt_end $b_tiltdir_end $b_origin_step]
   } elseif { $seqrtt eq "tr_t_r"} {
    set fnn [genbtrans $fnn $b_transstep $b_Ntrans]
    set fnn [genbt $fnn $b_Ntilt $b_tilt_end $b_tiltdir_end $b_origin_step]
    set fnn [genbr $fnn $b_stepr $b_Nrot]
   } else {
    error "Internal error: invalid sequence for animation"
   }  
   set headst 1 
 
  } else {
   incr fnn
  }

  if { $seqrtt eq "r_t_tr"} {
   set fnn [genhr $fnn $h_stepr $h_Nrot]
   set fnn [genht $fnn $h_Ntilt $h_tilt_end $h_tiltdir_end $h_origin_step]
   set fnn [genhtrans $fnn $h_transstep $h_Ntrans]
  } elseif { $seqrtt eq "t_r_tr"} {
   set fnn [genht $fnn $h_Ntilt $h_tilt_end $h_tiltdir_end $h_origin_step]
   set fnn [genhr $fnn $h_stepr $h_Nrot]
   set fnn [genhtrans $fnn $h_transstep $h_Ntrans]
  } elseif { $seqrtt eq "tr_r_t"} {
   set fnn [genhtrans $fnn $h_transstep $h_Ntrans]
   set fnn [genhr $fnn $h_stepr $h_Nrot]
   set fnn [genht $fnn $h_Ntilt $h_tilt_end $h_tiltdir_end $h_origin_step]
  } elseif { $seqrtt eq "tr_t_r"} {
   set fnn [genhtrans $fnn $h_transstep $h_Ntrans]
   set fnn [genht $fnn $h_Ntilt $h_tilt_end $h_tiltdir_end $h_origin_step]
   set fnn [genhr $fnn $h_stepr $h_Nrot]
  } else {
   error "Internal error: invalid sequence for animation"
  } 
 
  # start by rotating and tilting 
  lassign [::RADTOOL::get_vectors $MOLID $b_res0 $b_res1 $b_res2] b_V3 b_V1 b_V2 
  # these are the rotations based on the e coli structure
  # update for body and head.  this only has to be done one time to establish the coordinate system in the reference orientation. we don't need to update the head reference, since we will rotate the head first.
  set b_zero_tilt_dir [ ::RADTOOL::get_one_vector $MOLID $body_tilt0_0 $body_tilt0_1 ]
  set b_newone [vecnorm [veccross $b_V1 $b_zero_tilt_dir]]
  set b_newzero [vecnorm  [veccross $b_newone $b_V1]]
  
  # rotate about the axis of rotation.  
  set b_rot_vec $b_V1
 
  # start by rotating and tilting 
  lassign [::RADTOOL::get_vectors $MOLID $h_res0 $h_res1 $h_res2] h_V3 h_V1 h_V2 
  # these are the rotations based on the e coli structure
  # update for body and head.  this only has to be done one time to establish the coordinate system in the reference orientation. we don't need to update the head reference, since we will rotate the head first.
  set h_zero_tilt_dir [ ::RADTOOL::get_one_vector $MOLID $head_tilt0_0 $head_tilt0_1 ]
  set h_newone [vecnorm [veccross $h_V1 $h_zero_tilt_dir]]
  set h_newzero [vecnorm  [veccross $h_newone $h_V1]]
  
  # rotate about the axis of rotation.  
  set h_rot_vec $h_V1
 
  set endframe $fnn 
  if {$modn == 0 } {
   for {set ir 0} {$ir < $endframe} {incr ir} {
    set h_vec_end_rot [vecadd $rtd($ir,h_origin) $h_rot_vec]
    set b_vec_end_rot [vecadd $rtd($ir,b_origin) $b_rot_vec]

    animate dup frame 0 $MOLID
    $SELIDALL frame last
    $SELIDHEAD frame last
    moveelement $SELIDHEAD $rtd($ir,h_origin) $h_vec_end_rot $rtd($ir,h_rot) $rtd($ir,h_tiltd) $rtd($ir,h_tilt) $rtd($ir,h_trans) $h_newone $h_newzero  
    moveelement $SELIDALL $rtd($ir,b_origin) $b_vec_end_rot $rtd($ir,b_rot) $rtd($ir,b_tiltd) $rtd($ir,b_tilt) $rtd($ir,b_trans) $b_newone $b_newzero 
   }
   mol delete $MOLIDMOD 
  } elseif {$modn == 1} {
   for {set ir 0} {$ir < $endframe} {incr ir} {
    set h_vec_end_rot [vecadd $rtd($ir,h_origin) $h_rot_vec]
    set b_vec_end_rot [vecadd $rtd($ir,b_origin) $b_rot_vec]
    animate dup frame 0 $MOLIDMOD
    $SELIDALLMOD frame last
    $SELIDHEADMOD frame last 
    moveelement $SELIDHEADMOD $rtd($ir,h_origin) $h_vec_end_rot $rtd($ir,h_rot) $rtd($ir,h_tiltd) $rtd($ir,h_tilt) $rtd($ir,h_trans) $h_newone $h_newzero  
    moveelement $SELIDALLMOD $rtd($ir,b_origin) $b_vec_end_rot $rtd($ir,b_rot) $rtd($ir,b_tiltd) $rtd($ir,b_tilt) $rtd($ir,b_trans) $b_newone $b_newzero 
   }
   mol delete $MOLID 
  } else {
   error "internal error: invalid modn.  Please report issue to RADtool developers"
  }
  animate speed 0.9
  animate goto 0
  animate style Rock
  animate forward
  ::RADTOOL::orientribosome
 }

 # procs for generating sequences of angles

 proc initrtd {ONLYSSU b_rot_begin b_tilt_begin b_tiltdir_begin b_translate_begin b_origin_begin h_rot_begin h_tilt_begin h_tiltdir_begin h_translate_begin h_origin_begin} {
  variable rtd
  set fnn 0
  set rtd($fnn,h_rot)    $h_rot_begin 
  set rtd($fnn,h_tilt)   $h_tilt_begin
  set rtd($fnn,h_tiltd)  $h_tiltdir_begin
  set rtd($fnn,h_trans)  $h_translate_begin 
  set rtd($fnn,h_origin) $h_origin_begin

  # these are here to avoid an access error when only visualizing the head
  set rtd($fnn,b_rot)    0 
  set rtd($fnn,b_tilt)   0
  set rtd($fnn,b_tiltd)  0 
  set rtd($fnn,b_trans)  {0 0 0}
  set rtd($fnn,b_origin) {0 0 0} 

  if {$ONLYSSU == 1} { 
  # make a sequence of values for later 
   set rtd($fnn,b_rot)    $b_rot_begin
   set rtd($fnn,b_tilt)   $b_tilt_begin
   set rtd($fnn,b_tiltd)  $b_tiltdir_begin
   set rtd($fnn,b_trans)  $b_translate_begin
   set rtd($fnn,b_origin) $b_origin_begin
  } 

  return $fnn 
 }

 proc genbr  {fnn b_stepr b_Nrot} {
  variable rtd
  set pframe [expr $fnn-1]
  set br_begin $rtd($pframe,b_rot)
  if { $b_Nrot > 0 } {
   puts "    Body rotation begins with frame $fnn"
   for {set ir 1} {$ir < $b_Nrot+1} {incr ir} { 
    set rtd($fnn,b_rot)    [expr $br_begin+$b_stepr*$ir ]
    set rtd($fnn,b_tilt)   $rtd($pframe,b_tilt)  
    set rtd($fnn,b_tiltd)  $rtd($pframe,b_tiltd)
    set rtd($fnn,b_trans)  $rtd($pframe,b_trans)  
    set rtd($fnn,b_origin) $rtd($pframe,b_origin) 
                                               
    set rtd($fnn,h_rot)    $rtd($pframe,h_rot)      
    set rtd($fnn,h_tilt)   $rtd($pframe,h_tilt)   
    set rtd($fnn,h_tiltd)  $rtd($pframe,h_tiltd)  
    set rtd($fnn,h_trans)  $rtd($pframe,h_trans)    
    set rtd($fnn,h_origin) $rtd($pframe,h_origin) 
 
    incr fnn
   }
  }
  return $fnn
 }

 proc genbt {fnn b_Ntilt b_tilt_end b_tiltdir_end b_origin_step} {
  variable pi
  variable rtd
  variable deg2rad
  variable rad2deg
  set pframe [expr $fnn-1]
  set b_tilt_begin $rtd($pframe,b_tilt)
  set b_tiltdir_begin $rtd($pframe,b_tiltd)

  set pframe [expr $fnn-1]
  set bo_begin $rtd($pframe,b_origin)

  if {$b_Ntilt > 0} {
   puts "    Change in body tilt starts at frame $fnn"
   # set R axis for start
   set Rb  [list [expr cos($b_tiltdir_begin) ]  [expr sin($b_tiltdir_begin)] 0]
   set Rb [vecnorm [vecadd [vecscale [expr sin($b_tilt_begin*$deg2rad)] $Rb] [vecscale [expr cos($b_tilt_begin*$deg2rad)] {0 0 1}]]]

   set Re [list [ expr cos($b_tiltdir_end) ] [expr sin($b_tiltdir_end)] 0]
   set Re [vecnorm [vecadd [vecscale [expr sin($b_tilt_end*$deg2rad)] $Re] [vecscale [expr cos($b_tilt_end*$deg2rad)] {0 0 1}]]]
   set tiltstep [expr 1.0/$b_Ntilt]
   for {set ir 1} {$ir < $b_Ntilt+1} {incr ir} {
    set Rv [vecnorm [vecadd [vecscale [expr ($ir)*$tiltstep] $Re ]  [vecscale [expr (($b_Ntilt-$ir)*$tiltstep)] $Rb ] ] ]
    set tilt [expr acos([lindex $Rv 2])*$rad2deg]
    set tiltdir [veccross $Rv {0 0 1}]
    if {[veclength $tiltdir] > 0} { 
     set tiltdir [vecnorm $tiltdir]
     set phi [expr acos([lindex $tiltdir 0])*$rad2deg]
     if { [lindex $tiltdir 1] < 0 } {
      set phi [expr (-1)*$phi]
     }
     set phi [expr $deg2rad*($phi+90.0)]
    } else {
     set phi 0
    }
    set rtd($fnn,b_rot)    $rtd($pframe,b_rot)
    set rtd($fnn,b_tilt)   $tilt
    set rtd($fnn,b_tiltd)  $phi
    set rtd($fnn,b_trans)  $rtd($pframe,b_trans)
    set rtd($fnn,b_origin) [vecadd $bo_begin [vecscale $ir $b_origin_step]]
   
    set rtd($fnn,h_rot)    $rtd($pframe,h_rot)   
    set rtd($fnn,h_tilt)   $rtd($pframe,h_tilt)  
    set rtd($fnn,h_tiltd)  $rtd($pframe,h_tiltd) 
    set rtd($fnn,h_trans)  $rtd($pframe,h_trans)  
    set rtd($fnn,h_origin) $rtd($pframe,h_origin)
    incr fnn
   } 
  }
  return $fnn
 }

 proc genbtrans {fnn b_transstep b_Ntrans} {
  variable rtd
  set pframe [expr $fnn-1]
  set btrans_begin $rtd($pframe,b_trans)

  if {$b_Ntrans >0} {
   puts "    Body translation starts at frame $fnn"
   for {set ir 1} {$ir < $b_Ntrans+1} {incr ir} { 
    set rtd($fnn,b_rot)    $rtd($pframe,b_rot)  
    set rtd($fnn,b_tilt)   $rtd($pframe,b_tilt) 
    set rtd($fnn,b_tiltd)  $rtd($pframe,b_tiltd)
    set rtd($fnn,b_trans)  [vecadd $btrans_begin [vecscale $ir $b_transstep]]
    set rtd($fnn,b_origin) $rtd($pframe,b_origin)
  
    set rtd($fnn,h_rot)    $rtd($pframe,h_rot)   
    set rtd($fnn,h_tilt)   $rtd($pframe,h_tilt)  
    set rtd($fnn,h_tiltd)  $rtd($pframe,h_tiltd) 
    set rtd($fnn,h_trans)  $rtd($pframe,h_trans)  
    set rtd($fnn,h_origin) $rtd($pframe,h_origin)
    incr fnn
   }
  }
  return $fnn
 }

 proc genhr  {fnn h_stepr h_Nrot} {
  variable rtd
  variable headst
  set pframe [expr $fnn-1]
  set hr_begin $rtd($pframe,h_rot)

  if { $h_Nrot > 0 } {
   puts "    Head rotation begins with frame $fnn"
   for {set ir $headst} {$ir < $h_Nrot+1} {incr ir} { 
    set rtd($fnn,b_rot)    $rtd($pframe,b_rot)   
    set rtd($fnn,b_tilt)   $rtd($pframe,b_tilt)  
    set rtd($fnn,b_tiltd)  $rtd($pframe,b_tiltd) 
    set rtd($fnn,b_trans)  $rtd($pframe,b_trans) 
    set rtd($fnn,b_origin) $rtd($pframe,b_origin)
 
    set rtd($fnn,h_rot)    [expr $hr_begin+$h_stepr*$ir ]  
    set rtd($fnn,h_tilt)   $rtd($pframe,h_tilt)  
    set rtd($fnn,h_tiltd)  $rtd($pframe,h_tiltd)   
    set rtd($fnn,h_trans)  $rtd($pframe,h_trans)  
    set rtd($fnn,h_origin) $rtd($pframe,h_origin)
 
    incr fnn
   }
  }
  set headst 1
  return $fnn
 }

 proc genht {fnn h_Ntilt h_tilt_end h_tiltdir_end h_origin_step} {
  variable pi
  variable rtd
  variable headst
  variable deg2rad
  variable rad2deg
  set pframe [expr $fnn-1]
  set h_tilt_begin $rtd($pframe,h_tilt)
  set h_tiltdir_begin $rtd($pframe,h_tiltd)

  set pframe [expr $fnn-1]
  set ho_begin $rtd($pframe,h_origin)

  if {$h_Ntilt > 0} {
   puts "    Change in head tilt starts at frame $fnn"
   # set R axis for start
   set Rb  [list [expr cos($h_tiltdir_begin) ]  [expr sin($h_tiltdir_begin)] 0]
   set Rb [vecnorm [vecadd [vecscale [expr sin($h_tilt_begin*$deg2rad)] $Rb] [vecscale [expr cos($h_tilt_begin*$deg2rad)] {0 0 1}]]]

   set Re [list [ expr cos($h_tiltdir_end) ] [expr sin($h_tiltdir_end)] 0]
   set Re [vecnorm [vecadd [vecscale [expr sin($h_tilt_end*$deg2rad)] $Re] [vecscale [expr cos($h_tilt_end*$deg2rad)] {0 0 1}]]]
   set tiltstep [expr 1.0/$h_Ntilt]
   for {set ir $headst} {$ir < $h_Ntilt+1} {incr ir} {
    set Rv [vecnorm [vecadd [vecscale [expr ($ir)*$tiltstep] $Re ]  [vecscale [expr (($h_Ntilt-$ir)*$tiltstep)] $Rb ] ] ]
    set tilt [expr acos([lindex $Rv 2])*$rad2deg]
    set tiltdir [veccross $Rv {0 0 1}]
    if {[veclength $tiltdir] > 0} { 
     set tiltdir [vecnorm $tiltdir]
     set phi [expr acos([lindex $tiltdir 0])*$rad2deg]
     if { [lindex $tiltdir 1] < 0 } {
      set phi [expr (-1)*$phi]
     }
     set phi [expr $deg2rad*($phi+90.0)]
    } else {
     set phi 0
    }

    set rtd($fnn,b_rot)    $rtd($pframe,b_rot)   
    set rtd($fnn,b_tilt)   $rtd($pframe,b_tilt)  
    set rtd($fnn,b_tiltd)  $rtd($pframe,b_tiltd) 
    set rtd($fnn,b_trans)  $rtd($pframe,b_trans) 
    set rtd($fnn,b_origin) $rtd($pframe,b_origin)

    set rtd($fnn,h_rot)    $rtd($pframe,h_rot) 
    set rtd($fnn,h_tilt)   $tilt
    set rtd($fnn,h_tiltd)  $phi
    set rtd($fnn,h_trans)  $rtd($pframe,h_trans)  
    set rtd($fnn,h_origin) [vecadd $ho_begin [vecscale $ir $h_origin_step]]
    incr fnn
   } 
  }
  set headst 1
  return $fnn
 }

 proc genhtrans {fnn h_transstep h_Ntrans} {
  variable rtd
  variable headst
  set pframe [expr $fnn-1]
  set htrans_begin $rtd($pframe,h_trans)

  if {$h_Ntrans >0} {
   puts "    Head translation starts at frame $fnn"
   for {set ir $headst} {$ir < $h_Ntrans+1} {incr ir} { 
    set rtd($fnn,b_rot)    $rtd($pframe,b_rot)    
    set rtd($fnn,b_tilt)   $rtd($pframe,b_tilt)   
    set rtd($fnn,b_tiltd)  $rtd($pframe,b_tiltd)  
    set rtd($fnn,b_trans)  $rtd($pframe,b_trans)  
    set rtd($fnn,b_origin) $rtd($pframe,b_origin) 
 
    set rtd($fnn,h_rot)    $rtd($pframe,h_rot)   
    set rtd($fnn,h_tilt)   $rtd($pframe,h_tilt)  
    set rtd($fnn,h_tiltd)  $rtd($pframe,h_tiltd) 
    set rtd($fnn,h_trans)  [vecadd $htrans_begin [vecscale $ir $h_transstep]] 
    set rtd($fnn,h_origin) $rtd($pframe,h_origin)
    incr fnn
   }
  }
  set headst 1
  return $fnn
 } 

 # procs for generating step sizes
 proc getrotstep {rot_begin rot_end Nrot} {
  if {$Nrot < 0} {
   puts "Number of rotation values can not be negative."
   return
  } elseif {$Nrot > 0} {
   set stepr [expr 1.00*($rot_end-$rot_begin)/$Nrot]
  } else {
   if {$rot_end != $rot_begin} {
    error "If making 0 reps along rotation, rot_end and rot_begin must be equal "
   }
   set stepr 0
  }
  return $stepr
 }
 
 proc gettiltstep {tilt_begin tilt_end Ntilt} {
  if { ($tilt_end > 180) || ($tilt_end < 0)} {
   puts "Tilt angle is only defined for 0 to 180 degrees"
   return
  }
  if { ($tilt_begin > 180) || ($tilt_begin < 0)} {
   puts "Tilt angle is only defined for 0 to 180 degrees"
   return
  }
  if {$Ntilt < 0} {
   puts "Number of tilt values can not be negative."
   return
  } elseif {$Ntilt > 0} {
   set stept [expr 1.00*($tilt_end-$tilt_begin)/$Ntilt]
  } else {
   if {$tilt_end != $tilt_begin}  {
    error "If making 1 rep along tilt, and tilt is initially non-zero, tilt_end and tilt_begin must be equal"
   }
   set stept 0
  }
  return $stept
 }
 
 proc getoriginstep {tilt_begin tilt_end origin_begin origin_end Ntilt} {
  if {$tilt_begin < 0.1 && $tilt_end < 0.1} {
   set stepsize [list 0 0 0] 
  } else {
   if {$tilt_begin < 0.1} {
    set origin_begin $origin_end
   }
   if {$tilt_end < 0.1} {
    set origin_end $origin_begin
   }
   set stepsize [vecscale [expr 1.0/$Ntilt] [vecsub $origin_end $origin_begin]]
  }
  return [list $origin_begin $origin_end $stepsize]
 }
 
 proc gettiltdstep {tiltdir_begin tiltdir_end Ntiltd} {
  variable pi
 
  if {$Ntiltd < 0} {
   puts "Number of tilt direction values can not be negative."
   return
  } elseif {$Ntiltd > 0} {
   while {$tiltdir_end-$tiltdir_begin > 180} {
    set tiltdir_begin [expr $tiltdir_begin+360]
   }
   while {$tiltdir_end-$tiltdir_begin < -180} {
    set tiltdir_end [expr $tiltdir_end+360]
   }
   set steptd [expr 1.00*($tiltdir_end-$tiltdir_begin)/$Ntiltd]
  } else {
   if {($tiltdir_end != $tiltdir_begin) && ($tilt_begin != 0)} {
    error "If making 1 rep along tilt direction, and initial tilt is non-zero, tilt_end and tilt_begin must be equal"
   }
   set steptd 0
  }
  set tiltdir_begin [expr $pi*$tiltdir_begin/180.0 ]
  set tiltdir_end [expr $pi*$tiltdir_end/180.0 ]
  set steptd [expr $pi*$steptd/180.0 ]
  return [list $tiltdir_begin $tiltdir_end $steptd]
 }
 
 proc gettransstep {translate_begin translate_end Ntrans} {
  if {$Ntrans < 0} {
   puts "Number of translation points can not be negative."
   return
  } elseif {$Ntrans > 0} {
   set transstep [vecscale [expr 1.0/$Ntrans] [vecsub $translate_end $translate_begin] ]
  } else {
   if {[veclength $translate] > 0} {
    error "If showing 0 reps along translation, the trans vector must be 0"
   }
  }
  return $transstep
 }
 
 proc moveelement {SELID origin vec_end_rot rotby tiltdir tiltby transby onevec zerovec} {
 
  # rotate
  $SELID move [trans bond $origin $vec_end_rot $rotby deg]
  # set tilt
  # determine where the line of nodes should be.  We will calculate this as a rotation about V1.  Thus, we will take sin(theta)V2+cos(theta)V3
  set lofnodes [vecadd [ vecscale [expr sin($tiltdir)] $onevec ]  [ vecscale [expr cos($tiltdir)] $zerovec ]]
  #tilt
  set rot_vec $lofnodes
  set vec_end [vecadd $origin $rot_vec]
  $SELID move [trans bond $origin $vec_end $tiltby deg]
  $SELID moveby $transby
 }

 proc checkfinalorientation {SSUONLY rigidbody rigidhead animation hs he} {
  animate goto end
  set selr [atomselect $rigidhead "resid $hs to $he"]
  set selan [atomselect $animation "resid $hs to $he"]
  set RMSD [measure rmsd $selan $selr]
  if {$RMSD > 0.0001} {
   puts stderr "When generating the animation, the final interpolated position of the head is not exactly aligned with the rigid-body approximation. It is possible that the RADtool script was corrupted, or there was a conflict with another library. Things may look a bit off. The RMSD was $RMSD"
  }
  $selan delete
  $selr delete
  if {! $SSUONLY } {
   set selr [atomselect $rigidbody "not resid $hs to $he"]
   set selan [atomselect $animation "not resid $hs to $he"]
   set RMSD [measure rmsd $selan $selr]
   if {$RMSD > 0.0001} {
    puts stderr "When generating the animation, the final interpolated position of the body is not exactly aligned with the rigid-body approximation. It is possible that the RADtool script was corrupted, or there was a conflict with another library. Things may look a bit off. The RMSD was $RMSD"
   }
   $selan delete
   $selr delete
  }
 }
 
 proc ERfromEuler {rot_begin tilt_begin tiltdir_begin rot_end tilt_end tiltdir_end} {
  lassign [genvecs $rot_begin $tilt_begin $tiltdir_begin] v1 v2 v3
  lassign [genvecs $rot_end $tilt_end $tiltdir_end] v4 v5 v6
  return [lindex [::RADTOOL::find_ER_between_2_systems $v1 $v2 $v3 $v4 $v5 $v6] 0]
 }
 
 proc genvecs {rot tilt tiltdir} {

  set vx {1 0 0}
  set vy {0 1 0}
  set vz {0 0 1}

  set onevec {1 0 0}
  set zerovec {0 1 0}

  set M [trans bond {0 0 0} {0 0 1} $rot deg]
  set vx [vectrans $M $vx]
  set vy [vectrans $M $vy]

  set lofnodes [vecadd [ vecscale [expr sin($tiltdir)] $onevec ]  [ vecscale [expr cos($tiltdir)] $zerovec ]]
  set M [trans bond {0 0 0} $lofnodes $tilt deg]

  set vx [vectrans $M $vx]
  set vy [vectrans $M $vy]
  set vz [vectrans $M $vz]
  
  return [list $vx $vy $vz]
 }

} 
