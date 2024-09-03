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

##########################################################################
#                   Disclaimer and Acknowledgement                       #   
# This script was gradually built over years, where various snippets     #
# have certainly been found through online references and then modified  #
# for our purposes. We want to thank all of the VMD users and developers #
# who have made so many nice suggestions available online!!              #
# We hope that this script may also be as useful to other users as they  #
# build elaborate Tcl scripts for VMD.                                   #
##########################################################################

namespace eval RADTOOL {
 namespace export RADtool 
 
 # Determine where the script is located, so the procs can find all 
 # necessary supporting files
 if {[info exists radenv(GUION)]} {
  puts "Loading the RADtool plugin"
  if { [vmdinfo version] ne "1.9.3" } {
   if { [vmdinfo version] eq "1.9.4" } {
    puts stderr "WARNING: YOU ARE USING VMD 1.9.4. THIS VERSION OF 
  VMD WAS NOT AVAILABLE AT THE TIME RADTOOL WAS RELEASED. SO, THERE 
  IS NO WAY TO KNOW IF THE PERFORMANCE WILL BE RELIABLE. CHECK AND 
  SEE IF A NEW VERSION OF RADTOOL IS AVAILABLE, OR PERHAPS USE VMD 1.9.3."
   } else {
    puts stderr "WARNING: YOU ARE USING VMD [vmdinfo version]!!!!
  RADtool HAS ONLY BEEN TESTED WITH VMD v1.9.3. WHILE WE HAVE 
  FOUND RADTOOL OFTEN WORKS WITH 1.9.4aX, PERFORMANCE WITH 
  PRE-RELEASE VERSIONS OF VMD IS LIKELY TO BE LESS STABLE.
  
  IF YOU RUN \"RADtool -test\", AND ALL TESTS ARE PASSED, THEN
  MOST FEATURES MAY BE WORKING PROPERLY. HOWEVER, THERE ARE 
  KNOWN ISSUES WITH ANY VERSION OTHER THAN 1.9.3. SEE README."
   }
  }
 }
 variable pi
 set pi 3.14159265358979323846
 variable radversion
 variable radenv
 variable RADREF
 variable NGUYENREF
 variable labels
 variable fullfilenames
 variable RTD 
 set radenv(RADtoolhistory) {}

 #################### MAIN PROC #####################
 
 proc RADtool args {
 
  variable pi 
  variable radenv
  variable RTD 

  if { ! [info exists radenv(ROTTEST)] } {
   splashmessage "stdout"
  }
 
  # Since we only consider rRNA, let's make a list of stuff to ignore
  set radenv(NOTSTUFF) "not water and not ions and not amino" 
  #process arguments/options
  if {[catch {
   array set options [parse_args $args]
  } errorMessage] != 0} {
   puts stderr "Error:  $errorMessage\n\nExiting without completing.\n"
   if {[info exists radenv(GUION)]} {
    return 1
   } else {
    return 
   }
  } elseif {[info exists options(HELP)]} {
   # already printed the help info.  This just ensures that everything quits without proceeding
  } elseif {[info exists options(MOREINFO)]} {
   description
  } elseif {[info exists options(TESTON)] || [info exists options(TESTALLON)] || [info exists options(TESTTRAJON)] } {
   if {[info exists options(TESTALLON)] } {
    ::RTTESTS::RT-tests all
   } elseif {[info exists options(TESTON)] } {
    ::RTTESTS::RT-tests
   } elseif {[info exists options(TESTTRAJON)] } {
    ::RTTESTS::RT-tests trajonly
   }
  } elseif {[info exists options(ANIMATE2)] } {
   ::RTANIMATE::RT-animate2 
  } elseif {[info exists options(DOWNLOAD)] || [info exists options(BUNDLE)]} {
   set val [::RTBUNDLE::RT-bundle $args]
   if { [info exists radenv(ROTTEST)] } {
    return $val
   }
   puts "\nCompleted analysis of all ribosomes identified in the bundle."
  } else {
   set val [RT-master options]
   return $val
   
  } 
 }
 
 ##############END OF MAIN PROC #####################
 proc AddToRTD {args} {
  # adds data to the RadToolData (RTD) array
  # args is just a list of variable names
  # upvar retrieves the values of the variables and then RTD is updated
  variable RTD
  foreach val $args {
   upvar $val theval
   set RTD($val) $theval
  }
 } 

 proc RT-master {optin} {
 
  upvar $optin options
  variable radenv
  variable radversion
  variable RTD
 
  # set reference structural model names
  array set options [setreferences options] 
 
  # open all relevant file handles for writing.
  if {[catch {
   array set options [openoutputs options] 
  } errorMessage] != 0} {
   puts stderr "Error: $errorMessage\n\nExiting without completing.\n"
   if {[info exists radenv(GUION)]} {
    return 1
   } else {
    return 
   }
  }
  
  # Before loading anything, let's do our alignments, if selected.
  # Aligning earlier prevents alignment steps from getting bogged down by 
  # lots of loaded sequences.
  set RTD(nchangeSrRNA) "N/A"
  set RTD(nchangeLrRNA) "N/A"
  set RTD(nalignSrRNA) "N/A"
  set RTD(nalignLrRNA) "N/A"
  if {[info exists options(STAMP)]  || [info exists options(READALIGNMENT)] } {
   # first, perform a single stamp alignment, we don't want to perturb 
   # the input model, so we will simply move our references.
   if {[catch {
    set alignmaps [align_sequences options]
    AddToRTD alignmaps 
    } errorMessage] != 0} {
    if { ![info exists options(HEADLAST)]  } {
     set FHmess "\nPerhaps you should try to specify the head domain residues with -h_f and -h_l\n"
    } else {
     set FHmess ""
    }
    if {$options(ERROR) ne "stderr"} {
     puts $options(ERROR) "\n\nALIGNMENT ERROR: $errorMessage\n$FHmess\nExiting without completing.\n"
    }
    closelogs options
    error "\n\nALIGNMENT ERROR: $errorMessage\n$FHmess\nExiting without completing.\n"
   }
   set hf [lindex $alignmaps 2] 
   set hl [lindex $alignmaps 3]
   if { [llength $hf] != 0} {
    set options(HEADFIRST) $hf
   } 
   if { [llength $hl] != 0} {
    set options(HEADLAST) $hl
   } 
  }
 
  # load reference configs and core atoms
  if {[catch {
   loadrefdata options
  } errorMessage] != 0} {
   puts $options(ERROR) "Error: $errorMessage\n\nExiting without completing.\n"
   closelogs options
   if {[info exists radenv(GUION)]} {
    return 1
   } else {
    return 
   }
  }
 
  # set the residues define rigid-body axes
  defineresidues 
  # assign IDs and load structures if necessary
  loadmodel options
 
  # this is also done in the sequence alignment steps. But, we'll repeat it here, just in case we are not aligning.
  if {[catch {
   checkcontiguous options  
  } errorMessage] != 0} {
   puts stderr "Error: $errorMessage\n\nExiting without completing.\n"
   if {[info exists radenv(GUION)]} {
    return 1
   } else {
    return 
   }
  }
  # if we did an alignment earlier, renumber at this time.
  if { [info exists alignmaps] } {
   align_renum options
  }
 
  # determine the common atoms between the model and reference
  docommonstuff options 

  # prune the cores, if necessary
  if { [info exists options(PRUNE)]} {
   pruneall options
  }
 
  # save the cores to file, if needed
  if { [info exists options(CORESOUT)]} {
   puts $options(CORESOUT) "Cores Output file
Generated by RADtool v$radversion, on [clock format [clock second] -format "%D %T"], using machine [info hostname]"
   if { ! [info exists options(ONLYSSU)] } {
    puts $options(CORESOUT) "LSUCORE $RTD(COMMON_LrRNA)" 
   }
   if { ! [info exists options(ONLYLSU)] } {
    puts $options(CORESOUT) "BODYCORE $RTD(COMMON_BODY)" 
    puts $options(CORESOUT) "HEADCORE $RTD(COMMON_HEAD)" 
   }
   close $options(CORESOUT)
  }
 
  if { [info exists options(LOADTRAJ)]} { 
   if { ! [info exists options(ONLYSSU)] } {
    puts $options(OUTPUT) "# sequence matching: $RTD(NUM_COMMON_LrRNA) large rRNA residues, $RTD(NUM_COMMON_BODY) small rRNA body residues, $RTD(NUM_COMMON_HEAD) small rRNA head residues"
   } else {
    puts $options(OUTPUT) "# sequence matching: $RTD(NUM_COMMON_BODY) small rRNA body residues, $RTD(NUM_COMMON_HEAD) small rRNA head residues"
   }
  }
  if { ! [info exists options(ONLYSSU)] && $RTD(NUM_COMMON_LrRNA) ==0 } {
   puts $options(ERROR) "\nError: Must have at least one common residue in the large subunit. Found: $RTD(NUM_COMMON_LrRNA) \n"
   closelogs options
   return
  }
  if { ! [info exists options(ONLYLSU)] && ($RTD(NUM_COMMON_BODY) ==0 || $RTD(NUM_COMMON_HEAD) ==0 ) } {
   puts $options(ERROR) "\nError: Must have at least one common residue in the small head and small body. Found: $RTD(NUM_COMMON_HEAD), $RTD(NUM_COMMON_BODY) \n";
   closelogs options
   return
  }
  
  # set the representations for all of the rigid bodies
  drawrefs options
  ##############DEFINE ATOM SELECTIONS ##################
  if { ! [info exists options(ONLYSSU)] } {
   set RTD(p0_LrRNA) [atomselect $RTD(ID_REF_LrRNA) "resid $RTD(COMMON_LrRNA) and name P and not insertion A to Z"]
   set RTD(p1_LrRNA) [atomselect $RTD(ID_LARGE) "chain \"$RTD(CHAIN_ID_LrRNA)\" and resid $RTD(COMMON_LrRNA) and name P and not insertion A to Z and not altloc B to Z"]
   set RTD(large) [atomselect $RTD(ID_LARGE) all]
   if { ! [info exists options(ONLYLSU)] } {
    set RTD(body_c_0) [get_xyz $RTD(ID_REF_SrRNA_1) $RTD(body1)]  
   } else {
    set RTD(body_c_0) 0
   }
  } else {
   set RTD(p0_LrRNA) 0
   set RTD(p1_LrRNA) 0
   # weRTD( only need this defined, but we shouldn't use it.
   set RTD(bodytvecs) { {1 0 0} {0 1 0} {0 0 1}}
   set RTD(large) 0
   set RTD(body_c_0) 0
  }
  
  if { ! [info exists options(ONLYLSU)] } {
   set RTD(p1_body) [atomselect $RTD(ID_SMALL) "chain \"$RTD(CHAIN_ID_SrRNA)\" and resid $RTD(COMMON_BODY) and name P and not insertion A to Z and not altloc B to Z"]
   set RTD(p0_body) [atomselect $RTD(ID_REF_SrRNA_1) "resid $RTD(COMMON_BODY) and name P and not insertion A to Z"]
   set RTD(p1_head) [atomselect $RTD(ID_SMALL) "chain \"$RTD(CHAIN_ID_SrRNA)\" and resid $RTD(COMMON_HEAD) and name P and not insertion A to Z and not altloc B to Z"]
   set RTD(p0_head) [atomselect $RTD(ID_REF_SrRNA_2) "resid $RTD(COMMON_HEAD) and name P and not insertion A to Z"]
   set RTD(ref16_1) [atomselect $RTD(ID_REF_SrRNA_1) "all"]
   set RTD(ref16_2) [atomselect $RTD(ID_REF_SrRNA_2) "all"]
   set RTD(small) [atomselect $RTD(ID_SMALL) "all"]
   ###########END ATOM SELECTIONS#########################
   # some body reference values
   lassign [get_vectors $RTD(ID_REF_SrRNA_1) $RTD(body0) $RTD(body1) $RTD(body2)]  RTD(ref_bond_body) RTD(ref_norv_body) RTD(ref_perp_body)
   set RTD(zero_body) [ get_one_vector $RTD(ID_REF_SrRNA_1) $RTD(body_tilt0_0) $RTD(body_tilt0_1) ]
  } else {
   set RTD(p0_body) 0
   set RTD(p1_body) 0
   set RTD(p0_head) 0
   set RTD(p1_head) 0
   set RTD(ref16_1) 0
   set RTD(ref16_2) 0
   set RTD(small) 0
   # some body reference values
   set RTD(ref_bond_body) {0 0 0} 
   set RTD(ref_norv_body) {0 0 0}
   set RTD(ref_perp_body) {0 0 0}
   set RTD(zero_body) {0 0 0}
  }
 
 
  if { [info exists options(ONLYLSU)] } {
   # just run the initial step of the angle routines, since we are only aligning the LSU
   set RTD(LSURMSD) [do_all_angles options 0]
   set LSURMSD [format %.2f $RTD(LSURMSD)]
   puts $options(OUTPUT) "RMSD of LSU core between ref. and model: $RTD(LSURMSD)" 
  } else {
 
   if {! [info exists options(isTraj)] } {
    set framen 0
    # if not a trajectory, just analyze the loaded frame
    set out [do_all_angles options $framen]
 
    findbothtrios 
 
    if { [info exists radenv(ROTTEST)] } {
     close $options(OUTPUT)
     return $out
    }
 
   } elseif { $options(TRAJSEGLENGTH) != -1 } {
    # analyze in pieces of length $options(TRAJSEGLENGTH)
    # if trajectory, clear the frames before loading new ones.
 
    set seglength $options(TRAJSEGLENGTH)
    animate delete  beg 0 end -1 skip 0 $RTD(ID_LARGE)
 
    for {set i 1} {$i < [llength $options(LOADTRAJ)]} {incr i} {
     # for each file, load in pieces
     set loadnumb 0
     set loadnume [expr $seglength-1]
     set nf $seglength
     while {$seglength == $nf} {
       # if true, then the last load was able obtain a full segment, so continue
 
      set tfile [lindex $options(LOADTRAJ) $i]
      puts $options(OUTPUT) "Loading $tfile (frames $loadnumb to $loadnume)..." 
      # clear frames 
      animate delete  beg 0 end -1 skip 0 $RTD(ID_LARGE)
      # load on traj file
      mol addfile $tfile first $loadnumb last $loadnume step $options(TRAJSTEP) waitfor -1 $RTD(ID_LARGE)
 
      set nf [molinfo $RTD(ID_LARGE) get numframes]
      # LOOP OVER TRAJECTORY
      for {set framen 0} {$framen < $nf} {incr framen} { 
       # update the frame
       $RTD(p1_LrRNA) frame $framen
       $RTD(p1_body) frame $framen
       $RTD(p1_head) frame $framen
 
       do_all_angles options $framen
 
      #end of trajectory loop
      }
 
     set loadnumb [expr $loadnumb+$seglength]
     set loadnume [expr $loadnume+$seglength]
 
     }
    }
    close $options(TRAJOUT)
 
   } else {
     
    # if trajectory, clear the frames before loading new ones.
    animate delete  beg 0 end -1 skip 0 $RTD(ID_LARGE)
 
    for {set i 1} {$i < [llength $options(LOADTRAJ)]} {incr i} {
 
     set tfile [lindex $options(LOADTRAJ) $i]
     puts $options(OUTPUT) "Loading $tfile..." 
     # clear frames 
     animate delete  beg 0 end -1 skip 0 $RTD(ID_LARGE)
     # load on traj file
     mol addfile $tfile first $options(TRAJFIRST) last $options(TRAJLAST) step $options(TRAJSTEP) waitfor -1 $RTD(ID_LARGE)
     set nf [molinfo $RTD(ID_LARGE) get numframes]
     # LOOP OVER TRAJECTORY
     for {set framen 0} {$framen < $nf} {incr framen} { 
      # update the frame
      $RTD(p1_LrRNA) frame $framen
      $RTD(p1_body) frame $framen
      $RTD(p1_head) frame $framen
 
      do_all_angles options $framen
 
     #end of trajectory loop
     }
    }
    close $options(TRAJOUT)
   }
  }
  if { ! [info exists radenv(ROTTESTTRAJ)] } {
   puts "\nRADtool completed analysis of model.\n\n"
  } 
  if {$options(OUTPUT) ne "stdout"} {
   puts $options(OUTPUT) "\nRADtool completed analysis of model.\n"
  }
  if {! [info exists options(LOADTRAJ)]} {
   finalwrapup options 
  }
  closelogs options
 
  # clean up selections
 
  if { ! [info exists options(ONLYSSU)] } {
   $RTD(p0_LrRNA) delete 
   $RTD(p1_LrRNA)  delete
   $RTD(large)  delete
  }
  
  if { ! [info exists options(ONLYLSU)] } {
   $RTD(p1_body)  delete
   $RTD(p0_body)  delete
   $RTD(p1_head)  delete
   $RTD(p0_head)  delete
   $RTD(ref16_1)  delete
   $RTD(ref16_2)  delete
   $RTD(small)  delete
  }
 
 }
 
 ############## END RT-MASTER###########################3
 
 proc splashmessage {fh} {
  variable radversion
  variable RADREF
  puts $fh "
 
       [clock format [clock second] -format "%D %T"] 
       hostname: [info hostname]

       Ribosome Angle Decomposition (RAD) Tool, Version $radversion

       Copyright (c) 2021 The RAD development team at
       The Center for Theoretical Biological Physics,
       Northeastern University

       This tool is written for use with VMD. It takes structure 
       files (PDB, or CIF) containing the large subunit (LSU) 
       and/or small subunit (SSU) as input and returns the rotation, 
       tilt, tilt directon and translation of the body and head 
       of the SSU structure, where angles are calculated as described in:

$RADREF

       For more information about this tool, see:
            http://radtool.org

       Please send feedback to:
            Paul C. Whitford
            Center for Theoretical Biological Physics
            Northeastern University
            p.whitford@northeastern.edu
"
 
 }
 
 proc parse_args {in} {
  variable radenv
  set in [join $in]
  # give conversion between flags and internal variables
  # For historical reasons, this routine is unnecessarily convoluted.  
  # It will be cleaned up later.
  set flags(l) "LARGEPDB"
  set flags(s) "SMALLPDB"
  set flags(lc) "CHAIN_LrRNA"
  set flags(sc) "CHAIN_SrRNA"
  set flags(savepdb) "SAVENAME"
  set flags(traj) "LOADTRAJ"
  set flags(stamp) "STAMP"
  set flags(align_out) "ALIGNOUT"
  set flags(cores_out) "CORESOUT"
  set flags(test) "TESTON"
  set flags(testall) "TESTALLON"
  set flags(testtraj) "TESTTRAJON"
  set flags(ssuonly) "ONLYSSU"
  set flags(lsuonly) "ONLYLSU"
  set flags(notall) "NOINCLUDEALL"
  set flags(h_f) "HEADFIRST"
  set flags(h_l) "HEADLAST"
  set flags(noprune) "NOPRUNE"
  set flags(nofindhead) "NOFINDHEAD"
  set flags(findhead2) "FINDHEAD2"
  set flags(pruneby) "PRUNEBY"
  set flags(o) "OUTPUT"
  set flags(e) "ERROR"
  set flags(t_first) "TRAJFIRST"
  set flags(t_last) "TRAJLAST"
  set flags(t_step) "TRAJSTEP"
  set flags(ot) "TRAJOUT"
  set flags(align_in) "READALIGNMENT"
  set flags(cores_in) "READCORES"
  set flags(t_seglength) "TRAJSEGLENGTH"
  set flags(dump) "DUMPDATA"
  set flags(overwrite) "OVERWRITE"
  set flags(animate) "ANIMATE"
  set flags(animate2) "ANIMATE2"
  set flags(precision) "PRECISION"
  set flags(bundle) "BUNDLE"
  set flags(download) "DOWNLOAD"
 
  # set up reverse lookup.
  foreach name [array names flags] {
     set flagname($flags($name)) "-$name" 
  }
 
  if {[llength $in] == 0} {
   printusage
   set called(HELP)  "0"
   return [array get called]
  } else {
   set vars {}
   set flagarray ""
   set argl 0
   set seenany 0
   # if -help given, ignore all other flags and just give a message
   foreach opt $in {
    if { [regexp -nocase "^-h(elp)?$" $opt ] == 1  } {
     printusage
     set called(HELP)  "0"
     return [array get called]
    }
    if { [regexp -nocase "^-moreinfo$" $opt ] == 1  } {
     set called(MOREINFO)  "0"
     return [array get called]
    }
   }
 
   foreach opt $in {
    set c [regexp -all "^-" $opt]
    if {$c > 0} {
     set seenany 1
     regsub {^-}  $opt {} opt
     set opt [string tolower $opt]
     if {![info exists flags($opt)]} {
      error "Option $opt is not supported.  Try -help to see all supported options."
     } elseif { [info exists called($flags($opt))]} {
      error "-$opt given twice"
     } else {
      set opt $flags($opt)
      set called($opt) "0"
     }
     if {$argl == 0} {
      set currop $opt
     } else {
      set called($currop)  $vars
      set currop $opt
     }
     set vars {}
    } else {
     if {$seenany == 0} {
      error  "Bare word given: $opt
 Use -help to see supported options"
     }
     lappend vars $opt
    }
   incr argl
   }
   set called($currop)  $vars
  }
 
  foreach singleflags {TESTON TESTALLON TESTTRAJON ANIMATE2} {
   if { [ info exists called($singleflags)] } {
    # if it is a single flag, make sure it is the only value of called
    if { [llength [array names called]] != 1} {
     error "$flagname($singleflags) can not be given with any other options"
    }
    return [array get called] 
   }
  }
  # check that single file names were provided.
  foreach name {LARGEPDB TRAJOUT SMALLPDB READALIGNMENT READCORES CORESOUT} {
   if {[info exists called($name)]  && [llength $called($name)] != 1} {
    error  "Must provide a single file with option $flagname($name)"
   }
  }
 
  # check that single values were provided.
  foreach name {HEADFIRST HEADLAST TRAJSEGLENGTH CHAIN_LrRNA CHAIN_SrRNA PRECISION ANIMATE} {
   if {[info exists called($name)]  && [llength $called($name)] != 1} {
    error  "Must provide a single value with option $flagname($name)"
   }
  }

  foreach name [array names called] {
   set called($name) [join $called($name)]
  }

  if { [info exists called(READCORES)] } {
   set called(NOINCLUDEALL) 0
   puts "\nNote: $flagname(READCORES) given. Automatically turning on $flagname(NOINCLUDEALL) and $flagname(NOPRUNE)\n"
   set called(NOPRUNE) 0
  }
 
  # unless -hf and -hl are given with -stamp, then turn on find head
  if { [info exists called(STAMP)] && ![info exists called(HEADFIRST)] && ![info exists called(HEADLAST)] && ![info exists called(NOFINDHEAD)]} {
   set called(FINDHEAD) 0
  }
 
  # PRUNE and ALL are on, by default 
  if { ! [info exists called(NOPRUNE)] } {
   set called(PRUNE) 0
  }
  if { [info exists called(PRUNEBY)] } {
   if { ! [ string is double $called(PRUNEBY)] && ! [ string is integer $called(PRUNEBY)] } {
    error "Value given with $flagname(PRUNEBY) must be an integer or double.  Found \"$called(PRUNEBY)\""
   }
  } elseif { [info exists called(PRUNE)]} {
   set called(PRUNEBY) 2
  }
  if { ! [info exists called(NOINCLUDEALL)] } {
   set called(INCLUDEALL) 0
  }
  # set default precision
  if { ! [info exists called(PRECISION)] } {
   set called(PRECISION) 1
  }
 
  # check for required flags
  if { ! [info exists called(LOADTRAJ)] && ! [info exists called(BUNDLE)] && ! [info exists called(DOWNLOAD)] } {
   if { ! [info exists called(ONLYSSU)] } {
    foreach flag {LARGEPDB CHAIN_LrRNA } {
     if { ! [info exists called($flag)] } {
      error "Since you are not using -traj, -bundle, -download or -SSUonly, then $flagname($flag) is mandatory. To see description of options, use -h"
     }
    }
   }
   if { ! [info exists called(ONLYLSU)] } {
    foreach flag {SMALLPDB CHAIN_SrRNA} {
     if { ! [info exists called($flag)] } {
      error "Since you are not using -traj, -bundle, -download, or -LSUonly, then $flagname($flag) is mandatory. To see description of options, use -h"
     }
    }
   }
  }
 
  # check for incompatible options
  foreach {o1 o2} {SAVENAME BUNDLE SAVENAME DOWNLOAD SAVENAME LOADTRAJ NOPRUNE PRUNEBY LOADTRAJ ONLYSSU LOADTRAJ ONLYLSU BUNDLE CHAIN_LrRNA BUNDLE CHAIN_SrRNA BUNDLE LARGEPDB BUNDLE SMALLPDB DOWNLOAD LARGEPDB DOWNLOAD SMALLPDB DOWNLOAD CHAIN_LrRNA DOWNLOAD CHAIN_SrRNA BUNDLE DOWNLOAD ONLYSSU LARGEPDB ONLYSSU CHAIN_LrRNA ONLYLSU SMALLPDB ONLYLSU CHAIN_SrRNA INCLUDEALL CORESIN LOADTRAJ LARGEPDB LOADTRAJ SMALLPDB READALIGNMENT STAMP LOADTRAJ STAMP LOADTRAJ ONLYSSU LOADTRAJ ONLYLSU ONLYLSU FINDHEAD2} {
   if {[info exists called($o1)] && [info exists called($o2)] } {
    # add specific messages before throwing the error.
    set moreinfo ""
    if {$o1 eq "STAMP" && $o2 eq "LOADTRAJ"} {
     set moreinfo "Please run STAMP using a single structure, and then read in the sequence mapping when analyzing a trajectory"
    } elseif {$o1 eq "LOADTRAJ" && ($o2 eq "LARGEPDB" || $o2 eq "SMALLPDB" ) } {
     set moreinfo "When loading a trajectory, provide the molecule and trajectory names using the following format: \n\n\t$flagname($o1) <structure (both subunits in one file)> <trajectory files>"
    }
    error "Can not use $flagname($o1) and $flagname($o2) together. $moreinfo"
   }
  }
 
  # check that options were not given to specific flags
  foreach name {STAMP DUMPDATA OVERWRITE } {
   if {[info exists called($name)] && [llength $called($name)] != 0 } {
    error   " $flagname($name) can not be given a value. Found $called($name)"
   } 
  }
 
  # check for dependent options
  foreach {o1 o2} {LOADTRAJ CHAIN_LrRNA LOADTRAJ CHAIN_SrRNA FINDHEAD STAMP TRAJFIRST LOADTRAJ TRAJLAST LOADTRAJ TRAJSTEP LOADTRAJ } {
   if { [info exists called($o1)] && ! [info exists called($o2)] } {
    error   "Can not use $flagname($o1) without $flagname($o2)"
   }
  }
 
  # IF NOT SET, USE SOME DEFAULTS
  foreach {o1 o2} {ERROR stderr OUTPUT stdout TRAJSEGLENGTH -1 TRAJFIRST 0 TRAJLAST -1 TRAJSTEP 1 TRAJOUT traj.angles.out} {
   if { ! [info exists called($o1)] } {
    set called($o1) $o2
   }
  }
 
  # misc checks
  if { [info exists called(HEADLAST)] && [info exists called(HEADFIRST)] && $called(HEADLAST) < $called(HEADFIRST) } {
   error   "first head residue ($called(HEADFIRST) must be less than the last head residue ($called(HEADLAST))"
  }
 
  if { [info exists called(LOADTRAJ)] } {
   set called(SMALLPDB) [lindex $called(LOADTRAJ) 0]
   set called(LARGEPDB) [lindex $called(LOADTRAJ) 0]
  } 
 
  if {[info exists called(SMALLPDB)] && [info exists called(LARGEPDB)] && [info exists called(CHAIN_LrRNA)] && [info exists called(CHAIN_SrRNA)] } {
   if {$called(CHAIN_LrRNA) eq $called(CHAIN_SrRNA) && $called(LARGEPDB) eq $called(SMALLPDB)} {
    error "Using a single structure file ($called(LARGEPDB)) with both subunit given the same chain ID ($called(CHAIN_LrRNA)) not supported. "  
   }
  }
 
  if {[info exists called(SAVENAME)]} {
   regsub -all {.pdb$}    $called(SAVENAME) "" called(SAVENAME)
  }
 
 # if -traj selected, do some checks
  if {[info exists called(LOADTRAJ)]} {
   set called(LOADTRAJ)  [cleanup_single $called(LOADTRAJ) ] 
   set pdbfile [lindex $called(LOADTRAJ) 0]
   if { [regexp -all "\.pdb$" $pdbfile] != 1 } {
    error "First file listed with -traj option ($pdbfile) does not have suffix \".pdb\"\n"
   }
   for {set i 0} {$i < [llength $called(LOADTRAJ)]} {incr i} {
    set tfile [lindex $called(LOADTRAJ) $i]
    if { ! [file exists $tfile] } {
     error "File \"$tfile\" not found."
    }
   }
  }
 
  if { [info exists called(LOADTRAJ)] } {
   if { [llength $called(LOADTRAJ)] > 2 && $called(TRAJFIRST) != 0 } {
    error "Can not use -traj with multiple files (more than just a structure file and one trajectory) and also use -t_first"
   }
   if { [llength $called(LOADTRAJ)] > 2 && $called(TRAJLAST) != -1 } {
    error "Can not use -traj with multiple files (more than just a structure file and one trajectory) and also use -t_last"
   }
  }
 
  if { [info exists called(CHAIN_LrRNA) ] } {
   if { [info exists called(LOADTRAJ)] } {
    set ftmp [lindex $called(LOADTRAJ) 0]
   } else {
    set ftmp $called(LARGEPDB)
   }
   if { [regexp -nocase -all "\.pdb$" $ftmp] == 1 } {
    set frmt "PDB"
   } elseif { [regexp -nocase -all "\.cif$" $ftmp] == 1 } {
    set frmt "CIF"
   } else {
    error "Only PDB and CIF files are currently supported.  $ftmp given."
   }
   if { [string length $called(CHAIN_LrRNA)] > 2 && $frmt eq "CIF" } {
    error "LSU chain ID can not have more than two characters when using a cif file. Found \"$called(CHAIN_LrRNA)\"\n"
   }
   if { [string length $called(CHAIN_LrRNA)] != 1 && $frmt eq "PDB" } {
    error "LSU chain ID must be single characters when using a PDB file. Found \"$called(CHAIN_LrRNA)\"\n"
   }
  }
 
  if { [info exists called(CHAIN_SrRNA) ] } {
   if { [info exists called(LOADTRAJ)] } {
    set ftmp [lindex $called(LOADTRAJ) 0]
   } else {
    set ftmp $called(SMALLPDB)
   }
   if { [regexp -nocase -all "\.pdb$" $ftmp] == 1 } {
    set frmt "PDB"
   } elseif { [regexp -nocase -all "\.cif$" $ftmp] == 1 } {
    set frmt "CIF"
   } else {
    error "Only PDB and CIF files are currently supported.  $ftmp given."
   }
   if { [string length $called(CHAIN_SrRNA)] > 2 && $frmt eq "CIF" } {
    error "SSU chain ID can not have more than two characters when using a cif file. Found \"$called(CHAIN_SrRNA)\"\n"
   }
   if { [string length $called(CHAIN_SrRNA)] != 1 && $frmt eq "PDB" } {
    error "SSU chain ID must be single characters when using a PDB file. Found \"$called(CHAIN_SrRNA)\"\n"
   }
  }
 
  # we will use the output file name to define a file prefix for temporary STAMP files
  set tmp [lindex [split $called(OUTPUT) "/"] end]
  regsub -all \\.  "ribrottmp$tmp" "" called(TMPNAME)

  return [array get called]
 }

 proc setrefslabels {} {
  variable NGUYENREF
  variable RADREF
  variable labels

  set labels(arl) "aligned reference: large"
  set labels(aib) "aligned/idealized body"
  set labels(aih) "aligned/idealized head"
  set labels(rbp) "ref. body position"
  set labels(rhp) "ref. head position"
  set labels(ani) "animation"


  set labels(brar) "body rot. axis: ref."  
  set labels(bram) "body rot. axis: model"
  set labels(brpr) "body rot. plane: ref"
  set labels(brpm) "body rot. plane: model"
  set labels(bzta) "body zero-tilt axis" 
  set labels(bta)  "body tilt axis" 

  set labels(hrar) "head rot. axis: ref."  
  set labels(hram) "head rot. axis: model"
  set labels(hrpr) "head rot. plane: ref"
  set labels(hrpm) "head rot. plane: model"
  set labels(hzta) "head zero-tilt axis" 
  set labels(hta)  "head tilt axis" 

  set labels(erbra) "E-R body rot. axis"
  set labels(erhra) "E-R head rot. axis"

  set labels(tv1body) "body: trans v1"
  set labels(tv2body) "body: trans v2"
  set labels(tv3body) "body: trans v3"
  set labels(tv1head) "head: trans v1"
  set labels(tv2head) "head: trans v2"
  set labels(tv3head) "head: trans v3"

  set labels(btran) "body translation"
  set labels(htran) "head translation"

  set RADREF "       Hassan A, Byju S, Freitas FC, Roc C, Pender N, Nguyen K, Kimbrough EM,
       Mattingly J, Gonzalez RL Jr, Oliveira RJ, Dunham CM,  Whitford PC.
       Ratchet, swivel, tilt and roll: A complete description of subunit
       rotation in the ribosome. Nucleic Acids Research. 2022.
       DOI:10.1093/nar/gkac1211"
  set NGUYENREF "       Nguyen K, Whitford PC. 
       Steric interactions lead to collective head tilting 
       motion in the ribosome during mRNA-tRNA translocation.  
       Nature Communications. 7: 10586, 2016.
       DOI: 10.1038/ncomms10538"

 }
 
 proc setreferences {optin} {
 
  # define our reference models
  variable radenv
  upvar $optin options
  set options(BODYREF) "4v9d.16S.body.pdb" 
  set options(HEADREF) "4v9d.16S.head.pdb" 
  set options(SMALLREF) "4v9d.16S.pdb" 
  set options(LARGEREF) "4v9d.23S.pdb" 
  set options(LSUCORE) "23S-core.resid" 
  set options(BODYCORE) "16S-core.body.resid" 
  set options(HEADCORE) "16S-core.head.resid" 
  set SCRIPTPATH $radenv(ROTATIONPATH)
  set options(SCRIPTPATH) "$SCRIPTPATH/../"
 
  # check that all required reference files are accessible.
  foreach file {LARGEREF SMALLREF BODYREF HEADREF} {
   if {! [file exists $options(SCRIPTPATH)/share/reference_models/$options($file)]} {
    error "Unable to locate reference file $options($file).  This typically means the environment is not set properly. See README for details."
   }
  }
 
  foreach file {LSUCORE BODYCORE HEADCORE } {
   if {! [file exists $options(SCRIPTPATH)/share/coredefs/$options($file)]} {
    error "Unable to locate core definition file $options($file).  This typically means the environment is not set properly. See README for details."
   }
  }
 
  return [array get options]
 }
 
 proc openoutputs {optin} {
  variable radversion 
  upvar $optin options
  set outstr $options(OUTPUT)
  if {$outstr ne "stdout"} {
   if { [file exists $outstr ] && ![info exists options(OVERWRITE) ]} {
    error "$outstr already exists."
   } else {
    set drname [file dirname $outstr]
    if {![ file writable $drname]} {
     set fulldir [file normalize $drname]
     error "Unable to open output file ($outstr) in directory $fulldir\nMust have write access to save to this directory"
    } 
    set options(OUTPUT) [open "$outstr" w]
    splashmessage $options(OUTPUT)
   }
  }
 
  set errstr $options(ERROR)
  if {$errstr ne "stderr"} {
   if {$errstr eq $outstr} {
    set options(ERROR) $options(OUTPUT)
   } elseif { [file exists $errstr ] && ![info exists options(OVERWRITE) ]} {
    error "$errstr already exists."
   } else {
    set drname [file dirname $errstr]
    if {![ file writable $drname]} {
     set fulldir [file normalize $drname]
     error "Unable to open stdout file ($errstr) in directory $fulldir\nMust have write access to save to this directory"
    } 
    set options(ERROR) [open "$errstr" w]
   }
  }
  
  if { [info exists options(CORESOUT)] } {
   if { [file exists $options(CORESOUT) ] && ![info exists options(OVERWRITE) ]} {
    error "$options(CORESOUT) already exists."
   } else {
    set drname [file dirname $options(CORESOUT)]
    if {![ file writable $drname]} {
     set fulldir [file normalize $drname]
     error "Unable to open cores output file ($options(CORESOUT)) in directory $fulldir\nMust have write access to save to this directory"
    } 
    set options(CORESOUT) [open "$options(CORESOUT)" w]
   }
  } 
 
  if { [info exists options(LOADTRAJ)] } {
   set options(isTraj) 1
   if {[file exists $options(TRAJOUT)] && ![info exists options(OVERWRITE) ] } {
    error "$options(TRAJOUT) already exists."
    #return [array get options]
   } else {
    set drname [file dirname $options(TRAJOUT)]
    if {![ file writable $drname]} {
     set fulldir [file normalize $drname]
     error "Unable to open trajectory output file ($options(TRAJOUT)) in directory $fulldir\nMust have write access to save to this directory"
    } 

    set options(TRAJOUT) [open "$options(TRAJOUT)" w]
   }
   puts $options(TRAJOUT) "# Angle Trajectory, generated by RADtool v$radversion, on [clock format [clock second] -format "%D %T"], using machine [info hostname]"
   if { [info exists options(DUMPDATA)] } {
    puts $options(TRAJOUT) "# body_rotation body_phi body_psi body_tilt body_tilt_direction  head_rotation head_phi head_psi head_tilt head_tilt_direction RMSD_LrRNA RMSD_body RMSD_head {body_tilt_vector} {head_tilt_vector} ERangle_body {ERvec_body} ERangle_head {ERvec_head} {{Mlarge}} {{Mbody}} {{Mhead}} "
   } else {
    puts $options(TRAJOUT) "# body_rotation body_tilt body_tilt_direction head_rotation head_tilt head_tilt_direction RMSD_LrRNA RMSD_body RMSD_head "
   }
  } 
 
  if {  [info exists options(STAMP)]  && ! [info exists options(READALIGNMENT)] } {
   if {[info exists options(ALIGNOUT)] } {
    if {[file exists $options(ALIGNOUT)]  && ![info exists options(OVERWRITE) ]  } {
     error  "$options(ALIGNOUT) already exists."
    } else {
     set drname [file dirname $options(ALIGNOUT)]
     if {![ file writable $drname]} {
      set fulldir [file normalize $drname]
      error "Unable to open alignment output file ($options(ALIGNOUT)) in directory $fulldir\nMust have write access to save to this directory"
     } 
     set options(ALIGNOUT) [open "$options(ALIGNOUT)" w]
    }
   }
  }
  return [array get options]
 }
 
 proc checkcontiguous {optin} {
 
  upvar $optin options
  variable RTD
  foreach i { ID_LARGE CHAIN_ID_LrRNA ID_SMALL CHAIN_ID_SrRNA } { 
   set $i $RTD($i)
  }

  if { ! [info exists options(ONLYSSU)] } {
   puts $options(OUTPUT) "Checking LSU for continuous atom/residue numbering."
   # first see if any non rRNA atoms have the same chain ID.
   # if they do, then set their chains to "null"
   findlargestblock options $ID_LARGE $CHAIN_ID_LrRNA 1
   # for the largest block, check if atoms and resids are continuous 
   checkcontiguousatomsress options $ID_LARGE $CHAIN_ID_LrRNA 
  }
 
  if { ! [info exists options(ONLYLSU)] } {
   puts $options(OUTPUT) "Checking SSU for continuous atom/residue numbering."
   findlargestblock options $ID_SMALL $CHAIN_ID_SrRNA 1 
   checkcontiguousatomsress options $ID_SMALL $CHAIN_ID_SrRNA 
  }
 }

 proc checkcontiguousatomsress {optin MOL_ID CHAIN_ID} {
  upvar $optin options
  variable radenv
  set SEL [atomselect $MOL_ID "chain \"$CHAIN_ID\" and $radenv(NOTSTUFF)"]
  set first [lindex [$SEL get serial] 0]
  set last [lindex [$SEL get serial] [expr [$SEL num]-1]]
  set NUM [$SEL num]
  # check that the selected block of atoms is contiguous.
  # since call findlargestblock previously, this should never be true
  # but, it is possible that someone will find a way to break earlier 
  # routines (or we make a change that bypasses it) so, we'll leave this 
  # check in place
  if {[$SEL num] != [expr $last-$first+1] && [$SEL num] > 0} {
   set lastname [lindex [$SEL get name] [expr [$SEL num]-1]]
   set firstname [lindex [$SEL get name] 0]
   set lastresname [lindex [$SEL get resname] [expr [$SEL num]-1]]
   set firstresname [lindex [$SEL get resname] 0]
   error "Internal error: Contact developers. Non-contiguous set of atoms found.  This should not happen at this point in the calculation."
  }
  $SEL delete
  # check that the residue number are ascending and have no jumps.
  set seqmessages {}
  set SEL [atomselect $MOL_ID "chain \"$CHAIN_ID\" and name P and $radenv(NOTSTUFF)"]
  set residlist [$SEL get resid]
  set resnamelist [$SEL get resname]
  set last [lindex $residlist 0]
  for {set i 1} { $i < [llength $residlist] } {incr i} {
   set cur [lindex $residlist $i]
   set diff [expr $cur-$last]
   if { $diff > 1 || $diff < 0} {
    set curname [lindex $resnamelist $i]
    set lastname [lindex $resnamelist [expr $i-1]]
    lappend seqmessages "$lastname$last-$curname$cur"
   }
   set last $cur
  }
  if { [llength $seqmessages] > 0 } {
   puts $options(OUTPUT) "\nThe following instances of non-sequential residue numbering were found:"
   foreach v $seqmessages {
    puts $options(OUTPUT) "    $v"
   }
    puts $options(OUTPUT) ""
  } 
 } 
 
 proc loadrefdata {optin} {
  variable labels

  upvar $optin options
  set SCRIPTPATH $options(SCRIPTPATH)
 
  if { ! [info exists options(ONLYSSU)] } {
   set ID_REF_LrRNA  [mol new "$SCRIPTPATH/share/reference_models/$options(LARGEREF)" waitfor -1] 
   mol rename $ID_REF_LrRNA $labels(arl)
  } else {
   set ID_REF_LrRNA -1
  }
 
  if { ! [info exists options(ONLYLSU)] } {
   if { [info exists options(BODYREF)] } {
    set ID_REF_SrRNA_1  [mol new "$SCRIPTPATH/share/reference_models/$options(BODYREF)" waitfor -1] 
    mol rename $ID_REF_SrRNA_1 $labels(aib)
   } else {
    set ID_REF_SrRNA_1 -1
   }
   set ID_REF_SrRNA_2  [mol new "$SCRIPTPATH/share/reference_models/$options(HEADREF)" waitfor -1] 
   mol rename $ID_REF_SrRNA_2 $labels(aih)
   if { ! [info exists options(ONLYSSU)] } {
   # if only measuring head, then the zero position of the body is not defined
    set ID_REF_SrRNA_3  [mol new "$SCRIPTPATH/share/reference_models/$options(BODYREF)" waitfor -1] 
    mol rename $ID_REF_SrRNA_3 $labels(rbp) 
   } else {
    set ID_REF_SrRNA_3 -1
   }
   set ID_REF_SrRNA_4  [mol new "$SCRIPTPATH/share/reference_models/$options(HEADREF)" waitfor -1] 
   mol rename $ID_REF_SrRNA_4 $labels(rhp) 
 
   set ID_REF_SrRNA_5  [mol new "$SCRIPTPATH/share/reference_models/$options(SMALLREF)" waitfor -1] 
   mol rename $ID_REF_SrRNA_5 $labels(ani)
  } else {
    set ID_REF_SrRNA_1 -1
    set ID_REF_SrRNA_2 -1
    set ID_REF_SrRNA_3 -1
    set ID_REF_SrRNA_4 -1
    set ID_REF_SrRNA_5 -1
  }
 
  if {[info exists options(READCORES)] } {
   # if we are supposed to read the core definitions, then read.
   puts $options(OUTPUT) "Reading core definitions from file $options(READCORES)"
   set coresin [open "$options(READCORES)" r]
   set CORE_LrRNA {}
   set CORE_body {}
   set CORE_head {}
   while { ! [eof $coresin]} {
    set line [gets $coresin]
    if {[lindex $line 0] eq "LSUCORE" } {
     if {[llength $CORE_LrRNA] != 0} {
      error "Multiple LSU rRNA definitions found in cores file"
     }
     set CORE_LrRNA [lrange $line 1 end]
     puts $options(OUTPUT) "    LSU: [llength $CORE_LrRNA] residues"
    } 
 
    if {[lindex $line 0] eq "BODYCORE" } {
     if {[llength $CORE_body] != 0} {
      error "Multiple SSU body rRNA definitions found in cores file"
     }
     set CORE_body [lrange $line 1 end]
     puts $options(OUTPUT) "    SSU body: [llength $CORE_body] residues"
    }
 
    if {[lindex $line 0] eq "HEADCORE" } {
     if {[llength $CORE_head] != 0} {
      error "Multiple SSU head rRNA definitions found in cores file"
     }
     set CORE_head [lrange $line 1 end]
     puts $options(OUTPUT) "    SSU head: [llength $CORE_head] residues"
     # since we always write the head last, we can terminate once it is found
    }
   }
   if {! [info exists options(ONLYLSU)] } {
    if {[llength $CORE_body] == 0} {
     error "When reading cores, did not find the SSU body"
    }
    if {[llength $CORE_head] == 0} {
     error "When reading cores, did not find the SSU head"
    }
   }  

   if { ! [info exists options(ONLYSSU)] } {
    if {[llength $CORE_LrRNA] == 0} {
     error "When reading cores, did not find the LSU head"
    }
   }
   close $coresin
  } else {
   # default to the Whitford et al 2013 core definitions
   set data [open "$SCRIPTPATH/share/coredefs/$options(LSUCORE)" r]
   set CORE_LrRNA [read $data]
   close $data
   
   set data [open "$SCRIPTPATH/share/coredefs/$options(BODYCORE)" r]
   set CORE_body [read $data]
   close $data
   
   set data [open "$SCRIPTPATH/share/coredefs/$options(HEADCORE)" r]
   set CORE_head [read $data]
   close $data
  }
 
  set NUM_LrRNA [llength $CORE_LrRNA]
  set NUM_BODY [llength $CORE_body]
  set NUM_HEAD [llength $CORE_head]
  AddToRTD ID_REF_LrRNA ID_REF_SrRNA_1 ID_REF_SrRNA_2   ID_REF_SrRNA_3 ID_REF_SrRNA_4  ID_REF_SrRNA_5 CORE_LrRNA NUM_LrRNA CORE_body NUM_BODY CORE_head NUM_HEAD 
 }
 
 proc do_all_angles {optin framen } {
 
  variable radenv
  variable labels
  variable RTD
  foreach i { ref16_1 ref16_2 p0_LrRNA p1_LrRNA p1_body p0_body p1_head p0_head zero_body large small ID_LARGE ID_SMALL CHAIN_ID_LrRNA CHAIN_ID_SrRNA NUM_LrRNA NUM_BODY NUM_HEAD NUM_COMMON_LrRNA NUM_COMMON_BODY NUM_COMMON_HEAD ID_REF_SrRNA_1 ID_REF_SrRNA_2 ID_REF_SrRNA_3 ID_REF_SrRNA_4  ID_REF_SrRNA_5 ref_norv_body ref_perp_body ref_bond_body body0 body1 body2 body_d body_c_0 head0 head1 head2 head_d head_tilt0_0 head_tilt0_1 nchangeSrRNA nchangeLrRNA nalignSrRNA nalignLrRNA absorig } { 
   set $i $RTD($i)
  }
  upvar $optin options
  set prec $options(PRECISION)
  set prec "%.${prec}f"
 
  # get the absolute origin for the body - use the rigid-body before fitting, since this is in the
  # reference frame of the LSU.  If we are only analyzing the head, then it doesn't matter if this 
  # is not the correct frame of reference
  if { ! [info exists options(ONLYLSU)] } {
   set body_absorig [get_xyz $ID_REF_SrRNA_1 $absorig]
  }
 # do angle calculations
 
  if { [info exists options(ONLYLSU)] } {
   puts $options(OUTPUT) "Will only align the LSU and quit."
   set LSURMSD [calc_body_angle options $p0_LrRNA $p1_LrRNA $p1_body $p0_body $zero_body $large $small $framen $ID_LARGE $ID_SMALL $ID_REF_SrRNA_1 $ID_REF_SrRNA_4 $ID_REF_SrRNA_5 $ref_norv_body $ref_perp_body $ref_bond_body $body0 $body1 $body2 $body_d $body_c_0 $ref16_1 $ref16_2 $head_tilt0_0 $head_tilt0_1 {0 0 0} ] 
   return $LSURMSD 
  }
  if { ! [info exists options(ONLYSSU)] } {
   lassign [calc_body_angle options $p0_LrRNA $p1_LrRNA $p1_body $p0_body $zero_body $large $small $framen $ID_LARGE $ID_SMALL $ID_REF_SrRNA_1 $ID_REF_SrRNA_4 $ID_REF_SrRNA_5 $ref_norv_body $ref_perp_body $ref_bond_body $body0 $body1 $body2 $body_d $body_c_0 $ref16_1 $ref16_2 $head_tilt0_0 $head_tilt0_1 {0 0 0} ] zero_head RMSD_LrRNA RMSD_body body_rotation body_phi body_psi body_tilt body_tilt_direction body_tilt_vector axes_body_ref axes_body_fit ERvec_body ERangle_body body_ref_ax_pos ref_norv_body body_fit_ax_pos fit_norv_body zero_body Mlarge Mbody deltaabs_body bodytvecs
  } else {
   lassign [calc_body_angle options $p0_LrRNA $p1_LrRNA $p1_body $p0_body $zero_body $large $small $framen $ID_LARGE $ID_SMALL $ID_REF_SrRNA_1 $ID_REF_SrRNA_4  $ID_REF_SrRNA_5 $ref_norv_body $ref_perp_body $ref_bond_body $body0 $body1 $body2 $body_d $body_c_0 $ref16_1 $ref16_2 $head_tilt0_0 $head_tilt0_1 $body_absorig] zero_head RMSD_body Mbody deltaabs_body
  }
  # get the new axes for our translations of the head.  Use the reference head position (i.e. aligned based on body)
  # use the rigid body of the body for the origin of the head
  set head_absorig [get_xyz $ID_REF_SrRNA_1 $absorig]
  lassign [calc_head_angle options $p1_head $p0_head $zero_head $ref16_2 $ID_REF_SrRNA_1 $ID_REF_SrRNA_2 $head0 $head1 $head2 $head_d $head_absorig] RMSD_head head_rotation head_phi head_psi head_tilt head_tilt_direction head_tilt_vector axes_head_ref axes_head_fit ERvec_head ERangle_head head_ref_ax_pos ref_norv_head head_fit_ax_pos fit_norv_head zero_head Mhead headreset deltaabs_head headtvecs 
 
  # end angle calculations
  # Write out results 
  if { ! [info exists options(ONLYSSU)] } {
   set RMSD_LrRNA [format %.2f $RMSD_LrRNA]
  } else {
   set RMSD_LrRNA 0
  }
 
  set RMSD_body [format %.2f $RMSD_body]
  set RMSD_head [format %.2f $RMSD_head]
 
  # print the results
  if { [info exists options(isTraj)] } {
   # if we are running a trajectory, then we must reset the rigid body idealized position.  If we are not, then we leave it for visualization purposes
   $ref16_2 move $headreset
   if { [info exists radenv(ROTTEST)] } {
    set prec "%.3f"
   } 
 
   # if a trajectory, give different output format
   if { [info exists options(DUMPDATA)] } {
    if { ! [info exists options(ONLYSSU)] } {
     puts $options(TRAJOUT) "$body_rotation $body_phi $body_psi $body_tilt $body_tilt_direction  $head_rotation $head_phi $head_psi $head_tilt $head_tilt_direction $RMSD_LrRNA $RMSD_body $RMSD_head {$body_tilt_vector} {$head_tilt_vector} $ERangle_body {$ERvec_body} $ERangle_head {$ERvec_head} {$Mlarge} {$Mbody} {$Mhead}"
    } else {
     puts $options(TRAJOUT) "$head_rotation $head_phi $head_psi $head_tilt $head_tilt_direction $RMSD_body $RMSD_head {$head_tilt_vector} $ERangle_head {$ERvec_head} {$Mbody} {$Mhead}"
    }
   } else {
    foreach name {body_rotation body_tilt_direction body_tilt head_rotation head_tilt_direction head_tilt ERangle_body ERangle_head} {
     set j [set "$name"]
     set $name [format $prec $j]
    }
    if { ! [info exists options(ONLYSSU)] } {
     puts $options(TRAJOUT) "$body_rotation $body_tilt $body_tilt_direction $head_rotation $head_tilt $head_tilt_direction $RMSD_LrRNA $RMSD_body $RMSD_head"
    } else {
     puts $options(TRAJOUT) "$head_rotation $head_tilt $head_tilt_direction $ERangle_head $RMSD_body $RMSD_head"
    }
   }
  } else {
   lassign [linenearest $head_ref_ax_pos $ref_norv_head  $head_fit_ax_pos $fit_norv_head] head_origin head_final
   if { [info exists options(ONLYSSU)] } {
    set body_origin {0 0 0}
    set body_final {0 0 0}
    set body_rotation 0
    set body_tilt 0 
    set body_tilt_direction 0
    set body_tilt_vector {1 1 1}
    set axes_body_ref {1 1 1}
    set axes_body_fit {2 2 2}
    set ERvec_body {2 2 2}
    set ERangle_body 1
    set body_phi 0
    set body_psi 0
   } else {
    lassign [linenearest $body_ref_ax_pos $ref_norv_body  $body_fit_ax_pos $fit_norv_body] body_origin body_final
   }
 
   if { ! [info exists options(isTraj)]  } {
 
    # BODY ARROWS 
    if { ! [info exists options(ONLYSSU)] } {
     if { $body_tilt >0.1 } {
      set axcb $body_origin
      set xdir $body_tilt_vector
     } else {
      set axcb $body_ref_ax_pos
      set xdir $zero_body
     }
     draw_fitted_axis $axcb $ref_norv_body 6 $labels(brar)  150 
     draw_fitted_axis $axcb $fit_norv_body 4 $labels(bram)  150
     draw_plane $axcb $ref_norv_body $xdir 6 $labels(brpr) 
     draw_plane $axcb $fit_norv_body $xdir 4 $labels(brpm) 
     if { $body_tilt > 0.1 } {
      draw_fitted_axis $body_origin $zero_body 6 $labels(bzta)  150
      draw_fitted_axis $body_origin $body_tilt_vector 7 $labels(bta) 150
     } else {
      puts $options(ERROR) "\nNote: Body tilt is less than 0.1 degree.  Will not try to visualize,\n\tsince the tilt direction has little meaning for small tilt values."
     }
    }
    # HEAD ARROWS 
    #only draw axes if not a trajectory
    if { $head_tilt >0.1 } {
     set axc $head_origin
     set xdir $head_tilt_vector
    } else {
     set axc $head_ref_ax_pos
     set xdir $zero_head
    }
    draw_fitted_axis $axc $ref_norv_head 6 $labels(hrar)
    draw_fitted_axis $axc $fit_norv_head 3 $labels(hram)
 
    draw_plane $axc $ref_norv_head $xdir 6 $labels(hrpr)
    draw_plane $axc $fit_norv_head $xdir 3 $labels(hrpm)
 
    if { $head_tilt >0.1 } {
     draw_fitted_axis $head_origin $zero_head 6 $labels(hzta)
     draw_fitted_axis $head_origin $head_tilt_vector 14 $labels(hta)
    } else {
     puts $options(ERROR) "\nNote: Head tilt angle is less than 0.1 degree, so the tilt axis is not well defined. Will not try to visualize,\n\tsince the tilt direction has little meaning for small tilt values."
    }
   } 
 
   ####### Find ER-associated translation############
   if { ! [info exists options(ONLYSSU)]} {
    if { $ERangle_body != 0 } {
     lassign [find_COR_with_minimum_absolute_translation $Mbody "BODY" [transidentity]] COR_body minimum_translation_body	
    } else {
     set COR_body "not defined"
     set minimum_translation_body {0 0 0}	
    }
    if { $ERangle_head != 0 } {
     lassign [find_COR_with_minimum_absolute_translation $Mhead "HEAD" $Mbody] COR_head minimum_translation_head
    } else {
     set COR_head "not defined"
     set minimum_translation_head {0 0 0}	
    }
   } else {
    set Mlarge "not defined"
    if { $ERangle_head != 0 } {
     lassign [find_COR_with_minimum_absolute_translation $Mhead "HEAD" [transidentity]] COR_head minimum_translation_head
    } else {
     set COR_head "not defined"
     set minimum_translation_head {0 0 0}	
    }
    set COR_body "not defined"
    set minimum_translation_body {0 0 0}	
   }
   ###################
   set body_translation [vecsub $body_final $body_origin]
   set head_translation [vecsub $head_final $head_origin]
   if {[info exists options(HEADFIRST)]} {
    set hf $options(HEADFIRST)
   } else {
    set hf 928
   }
   if {[info exists options(HEADLAST)]} {
    set hl $options(HEADLAST)
   } else {
    set hl 1389
   }
 
   if { ! [info exists options(ONLYSSU)]} {
    lappend radenv(RADtoolhistory) [list $body_rotation $body_tilt $body_tilt_direction $body_translation $body_origin $head_rotation $head_tilt $head_tilt_direction $head_translation $head_origin $hf $hl $ID_REF_SrRNA_1 $ID_REF_SrRNA_2 $CHAIN_ID_SrRNA $ID_LARGE $ID_SMALL]
   } else {
    lappend radenv(RADtoolhistory) [ list 0 0 0 {0 0 0} {0 0 0} $head_rotation $head_tilt $head_tilt_direction $head_translation $head_origin $hf $hl $ID_REF_SrRNA_1 $ID_REF_SrRNA_2 $CHAIN_ID_SrRNA -100 $ID_SMALL]
   }
 
   if { $COR_body ne "not defined" } {
    if { $ERangle_body < 1 } {
     puts $options(ERROR) "E-R angle for the body is small. While the direction is meaningful, the location of the axis can be noisy"
    }
    draw_fitted_axis $COR_body $ERvec_body "yellow" $labels(erbra)  125 
   }
   if { $COR_head ne "not defined" } {
    if { $ERangle_head < 1 } {
     puts $options(ERROR) "E-R angle for the head is small. While the direction is meaningful, the location of the axis can be noisy"
    }
    draw_fitted_axis $COR_head $ERvec_head "yellow" $labels(erhra) 125 
   }
 
   if { ! [info exists options(ONLYSSU)] } {
   ## convert to internal coordinate system
    if {[veclength $body_translation] > 0} {
     draw_fitted_axis $axcb $body_translation 12 $labels(btran) [veclength $body_translation] 10 
    } else {
     puts "NOTE: Body translation is 0, will not try to draw vectors"
    }
    set bt_or $body_translation 
    set body_translation [transvecstrans $body_translation $bodytvecs ]
    drawtransvecstrans $bodytvecs "body"
    set deltaabs_body [transvecstrans $deltaabs_body $bodytvecs ]
    set minimum_translation_body [transvecstrans $minimum_translation_body $headtvecs ]
   }
 
   if {[veclength $head_translation] > 0} {
    draw_fitted_axis $axc $head_translation 12 $labels(htran) [veclength $head_translation] 10
   } else {
     puts "NOTE: Head translation is 0, will not try to draw vectors"
   } 
   set ht_or $head_translation 
   set head_translation [transvecstrans $head_translation $headtvecs ]
   drawtransvecstrans $headtvecs "head"
   set deltaabs_head [transvecstrans $deltaabs_head $headtvecs ]
   set minimum_translation_head [transvecstrans $minimum_translation_head $headtvecs ]
 
   set out [print_results options $ID_LARGE $ID_SMALL $CHAIN_ID_LrRNA $CHAIN_ID_SrRNA $body_rotation $body_phi $body_psi $body_tilt $body_tilt_direction $body_tilt_vector $body_translation $body_origin $head_rotation $head_phi $head_psi $head_tilt $head_tilt_direction $head_tilt_vector $head_translation $head_origin $RMSD_LrRNA $RMSD_body $RMSD_head $NUM_LrRNA $NUM_COMMON_LrRNA $NUM_BODY $NUM_COMMON_BODY $NUM_HEAD $NUM_COMMON_HEAD $axes_body_ref $axes_body_fit $axes_head_ref $axes_head_fit $ERvec_body $ERangle_body $ERvec_head $ERangle_head $nalignLrRNA $nalignSrRNA $nchangeLrRNA $nchangeSrRNA $COR_body $minimum_translation_body $COR_head $minimum_translation_head $deltaabs_body $deltaabs_head $Mlarge $Mbody $Mhead]
   if {[info exists radenv(ROTTEST)]} {
    return $out
   }
  }
  # end of writing results
 
  if {[info exists options(ANIMATE)]} {
   # turn on animation
   ::RTANIMATE::RT-animate options $ID_REF_SrRNA_5 $ID_REF_SrRNA_1 $ID_REF_SrRNA_2 $options(ANIMATE) 20 20 10 20 20 10   
  }
  return
 }
 
 proc linenearest {p1 d1 p2 d2} {
  # takes equations of two lines: pi+t*di and returns points along each at which the lines are closest
 
  set dp [vecsub $p1 $p2]
  set dv [vecsub $d2 $d1]
  set dpsum [expr [lindex $dp 0]*[lindex $dv 0]]
  set dpsum [expr $dpsum + [lindex $dp 1]*[lindex $dv 1]]
  set dpsum [expr $dpsum + [lindex $dp 2]*[lindex $dv 2]]
  set dv2 [veclength2 [vecsub $d1 $d2]] 
  if {$dpsum == 0} {
   set t 0
  } else { 
   if {$dv2 == 0} {
    error "Trying to find intersection of two lines, but the lines are parallel. This is extremely unlikely to ever occur for a typical ribosome structure. Please report this issue to the Whitford Group."
   }
   set t [expr $dpsum/$dv2]
  }
  set c1 [vecadd $p1 [vecscale $t $d1]]
  set c2 [vecadd $p2 [vecscale $t $d2]]
  return [list $c1 $c2]
 }
 
 proc transdiff {orig p1 p2 angle rotvec} {
  # takes three points (origin, point before and point after), the rotation angle and axis.  Then calculates the associated translation when uisng that origin to rotate the point
  # in principle, we could use "origin" instead of center and offset.  But, for some reason, origin wasn't centering properly.  
  set t [trans center $orig offset $orig axis $rotvec $angle deg]
  set p1r [coordtrans $t $p1]
  set deltax [vecsub $p2 $p1r]
  return $deltax
 }
 
 proc align_sequences {optin} {
  # this either applies a sequence alignment, STAMP alignment, or reads in the alignments and renumbers residues.
 
  upvar $optin options
  set hf {}
  set hl {}
  if {[info exists options(READALIGNMENT)]} {
   ##################### READ ALIGNMENT ##############################
   puts $options(OUTPUT) "\nSequence mapping will be read from $options(READALIGNMENT)"
   set nalignSrRNA "N/A"
   set nalignLrRNA "N/A"
   # if we performed the alignment previously, we can just read the map.
   set mapin [open "$options(READALIGNMENT)" r]
   set foundlarge 0
   set foundsmall 0 
   set foundsmallh 0 
   while { ! [eof $mapin]} {
    set list [gets $mapin]
    if { [lindex $list 0] eq "large" &&  [lindex $list 1] eq "subunit:" } {
     if {$foundlarge == 1} {
      error "Multiple LSU definitions found in mapping file."
     }
     set foundlarge 1
     if { [info exists options(ONLYSSU)] } {
      puts $options(OUTPUT) "  Found LSU mapping in the align_in file, but -SSUonly was given.  Will ignore entry."
      continue
     } else {
      set list [lrange $list 2 end]
      if { [expr [llength $list] % 2] != 0 } {
       error "Alignment file has an odd number of residues numbers listed
 fo    the LSU. The file should list pairs of numbers (model and E. coli number)"
      }
      foreach {name val} $list {
       set LrRNAmap($name) $val
      }
      puts $options(OUTPUT) "    LSU: [array size LrRNAmap] residues found in mapping file"
     } 
     if { [info exists options(ONLYLSU)] } {
      #since we only need the LSU, we won't even look for the SSU definition
      break
     }
    }
 
    if {  [lindex $list 0] eq "small" &&  [lindex $list 1] eq "subunit:" } {
     if {$foundsmall == 1} {
      error "Multiple SSU definitions found in mapping file."
     }
     set foundsmall 1
     set list [lrange $list 2 end]
     if { [expr [llength $list] % 2] != 0 } {
      error "Alignment file has an odd number of residues numbers listed
 for the SSU. The file should list pairs of numbers (model and E. coli number)"
     }
     foreach {name val} $list {
      set SrRNAmap($name) $val
     }
     puts $options(OUTPUT) "    SSU: [array size SrRNAmap] residues found in mapping file\n"
    }

    if {  [lindex $list 0] eq "small" &&  [lindex $list 1] eq "head:" } {
     if {$foundsmallh == 1} {
      error "Multiple SSU head definitions found in mapping file."
     }
     set foundsmallh 1
     set list [lrange $list 2 end]
     if { [llength $list] != 3 } {
      error "Alignment file has incorrect format. Found \"small head:\" line, but it is not followed by \"N to M\" where N and M are the first and last head residue"
     }
     set hf [lindex $list 0]
     set hl [lindex $list 2]
     # we never read past this line.
    }
   }
   if { ! [info exists options(ONLYLSU)]} {
    if {$foundsmall != 1 } {
    error "Failed to find SSU mapping information in alignment file"
    }
    if {$foundsmallh != 1 } {
    error "Failed to find SSU head definition in alignment file"
    }
   }
   if {$foundlarge != 1 && ! [info exists options(ONLYSSU)]} {
    error "Failed to find LSU mapping information in alignment file"
   }
   close $mapin
  } elseif { [info exists options(STAMP)] }  {
   ################### PERFORM ALIGNMENT ############################
   # perform alignment
   set tlist [doseqstampalign options] 
   array set options [lindex $tlist 0]
   if { ! [info exists options(ONLYSSU)] } {
    array set LrRNAmap [lindex $tlist 1]
   } else {
    set LrRNAmap(0) 0
   }
   if { ! [info exists options(ONLYLSU)] } {
    array set SrRNAmap [lindex $tlist 2]
   } else {
    set SrRNAmap(0) 0
   }
  } else {
   error "Internal error: Invalid alignment approach called" 
  }
  return [list [array get LrRNAmap] [array get SrRNAmap] $hf $hl]
 }
 
 proc align_renum {optin } {
  # this either applies a sequence alignment, STAMP alignment, or reads in the alignments and renumbers residues.
  upvar $optin options
  variable RTD
  foreach i { alignmaps ID_LARGE CHAIN_ID_LrRNA ID_SMALL CHAIN_ID_SrRNA  } {
   set $i $RTD($i)
  }
  
  set nalignLrRNA 0
  set nchangeLrRNA 0
  if { ! [info exists options(ONLYSSU)] } {
   array set LrRNAmap [lindex $alignmaps 0]
   set nalignLrRNA [array size LrRNAmap]
   # go residue by residue, using residue (not resid) and reassign resid
   set selt [atomselect $ID_LARGE "chain \"$CHAIN_ID_LrRNA\" and name P" ]
   foreach res [$selt get index] {
    set rest [atomselect $ID_LARGE "index $res"] 
    set setl [$rest get resid] 
    set insl [$rest get insertion]
    lassign [cleanup "$setl $insl"] a b 
    if { [info exists LrRNAmap($a$b) ] } {
     $rest set resid $LrRNAmap($a$b)
     if { $LrRNAmap($a$b) != "$a$b"} {
      incr nchangeLrRNA
     }
    } else {
     $rest set resid -100
    } 
    $rest delete
   }
   $selt delete
  } 
 
  set nalignSrRNA 0
  set nchangeSrRNA 0
  if { ! [info exists options(ONLYLSU)] } {
   array set SrRNAmap [lindex $alignmaps 1]
   set nalignSrRNA [array size SrRNAmap]
   #map reference model numbering to the input model
   # how many resids are changed.
 
   # go residue by residue, using residue (not resid) and reassign resid
   set selt [atomselect $ID_SMALL "chain \"$CHAIN_ID_SrRNA\" and name P" ]
   foreach res [$selt get index] {
    set rest [atomselect $ID_SMALL "index $res"] 
    set setl [$rest get resid] 
    set insl [$rest get insertion]
    lassign [cleanup "$setl $insl"] a b 
    if { [info exists SrRNAmap($a$b) ] } {
     $rest set resid $SrRNAmap($a$b)
     if { $SrRNAmap($a$b) != "$a$b"} {
      incr nchangeSrRNA
     }
    } else {
     $rest set resid -100
    } 
    $rest delete
   }
   $selt delete
  }
  AddToRTD nalignLrRNA nalignSrRNA nchangeLrRNA nchangeSrRNA
 }
 
 proc docommonstuff {optin} {
 
  variable RTD
  foreach i { ID_REF_LrRNA ID_LARGE CHAIN_ID_LrRNA ID_REF_SrRNA_1 ID_REF_SrRNA_2 ID_SMALL CHAIN_ID_SrRNA CORE_LrRNA CORE_body CORE_head } { 
   set $i $RTD($i)
  }

  upvar $optin options
  # determine which residues are common to the model and reference
  if { [info exists options(INCLUDEALL)]} { 
   puts $options(OUTPUT) "Will consider all possible common residues for structural alignment"
   # here, we will just take all aligned residues, or just all residues
   if { ! [info exists options(ONLYSSU)] } {
    lassign [getallcommon "Large rRNA" $ID_REF_LrRNA $ID_LARGE $CHAIN_ID_LrRNA "1 to 3000" 1200 $options(OUTPUT)] NUM_COMMON_LrRNA COMMON_LrRNA  
   } else {
    set NUM_COMMON_LrRNA -1
    set COMMON_LrRNA -1
   }
   if { ! [info exists options(ONLYLSU)] } {
    lassign [getallcommon "Small rRNA body" $ID_REF_SrRNA_1 $ID_SMALL $CHAIN_ID_SrRNA "2 to 927 1390 to 1600" 400 $options(OUTPUT)] NUM_COMMON_BODY COMMON_BODY 
    lassign [getallcommon "Small rRNA head" $ID_REF_SrRNA_2 $ID_SMALL $CHAIN_ID_SrRNA "928 to 1389" 120 $options(OUTPUT)] NUM_COMMON_HEAD COMMON_HEAD 
   } else {
    set NUM_COMMON_BODY -1
    set COMMON_BODY -1
    set NUM_COMMON_HEAD -1
    set COMMON_HEAD -1
   }
  } else {
   puts $options(OUTPUT) "\nWill consider only pre-defined CORE residues for structural alignment"
   # since we are not saying to use all, then just use core atoms
   # find the core atoms that are present in the model
   if { ! [info exists options(ONLYSSU)] } {
    lassign [getcommon "Large rRNA" $ID_LARGE $CHAIN_ID_LrRNA $CORE_LrRNA 1200 $options(OUTPUT)] NUM_COMMON_LrRNA COMMON_LrRNA  
   } else {
    set NUM_COMMON_LrRNA -1
    set COMMON_LrRNA -1
   }
   if { ! [info exists options(ONLYLSU)] } {
    lassign [getcommon "Small rRNA body" $ID_SMALL $CHAIN_ID_SrRNA $CORE_body 400 $options(OUTPUT)] NUM_COMMON_BODY COMMON_BODY 
    lassign [getcommon "Small rRNA head" $ID_SMALL $CHAIN_ID_SrRNA $CORE_head 120 $options(OUTPUT)] NUM_COMMON_HEAD COMMON_HEAD 
   } else {
    set NUM_COMMON_BODY -1
    set COMMON_BODY -1
    set NUM_COMMON_HEAD -1
    set COMMON_HEAD -1
   }
  }
  if {$NUM_COMMON_LrRNA == 0} {
   error "Number of common LSU residues is 0. This typically means the input structure is not using E coli numbering. In those cases, it is necessary to use the -stamp option."
  }
  if {$NUM_COMMON_BODY == 0} {
   error "Number of common SSU body residues is 0. This typically means the input structure is not using E coli numbering. In those cases, it is necessary to use the -stamp option."
  }
  if {$NUM_COMMON_HEAD == 0} {
   error "Number of common SSU head residues is 0. This typically means the input structure is not using E coli numbering. In those cases, it is necessary to use the -stamp option."
  }
  AddToRTD NUM_COMMON_LrRNA COMMON_LrRNA NUM_COMMON_BODY COMMON_BODY NUM_COMMON_HEAD COMMON_HEAD
 }
 
 proc pruneall {optin} {
 
  upvar $optin options
  variable RTD
  foreach i { ID_LARGE CHAIN_ID_LrRNA COMMON_LrRNA ID_REF_LrRNA ID_SMALL CHAIN_ID_SrRNA COMMON_BODY COMMON_HEAD ID_REF_SrRNA_1 ID_REF_SrRNA_2 } { 
   set $i $RTD($i)
  }
  # iteratively align and remove atoms that differ after alignment, until all atoms are with 1A
  # give the molID, chain and atoms.  This way, we can copy and analyze them. Also give the ref structure ID, so we can align to it
  puts $options(OUTPUT) "pruning (cutoff $options(PRUNEBY) Angstroms) the list of residues that will be used for structure alignment."
  if { ! [info exists options(ONLYLSU)] } {
   puts $options(OUTPUT) "pruning the SSU head"
   lassign [prunecore options $ID_SMALL $CHAIN_ID_SrRNA $COMMON_HEAD $ID_REF_SrRNA_2] NUM_COMMON_HEAD COMMON_HEAD
   puts $options(OUTPUT) "done pruning the SSU head"
   puts $options(OUTPUT) "pruning the SSU body"
   lassign [prunecore options $ID_SMALL $CHAIN_ID_SrRNA $COMMON_BODY $ID_REF_SrRNA_1] NUM_COMMON_BODY COMMON_BODY
   puts $options(OUTPUT) "done pruning the SSU body"
  } else {
   set NUM_COMMON_BODY -1
   set COMMON_BODY -1
   set NUM_COMMON_HEAD -1
   set COMMON_HEAD -1
  }
  if { ! [info exists options(ONLYSSU)] } {
   puts $options(OUTPUT) "pruning the LSU"
   lassign [prunecore options $ID_LARGE $CHAIN_ID_LrRNA $COMMON_LrRNA $ID_REF_LrRNA] NUM_COMMON_LrRNA COMMON_LrRNA
   puts $options(OUTPUT) "done pruning the LSU"
  } else {
   set NUM_COMMON_LrRNA -1
   set COMMON_LrRNA -1
  }
  AddToRTD NUM_COMMON_LrRNA COMMON_LrRNA NUM_COMMON_BODY COMMON_BODY NUM_COMMON_HEAD COMMON_HEAD
 }
 
 proc drawrefs {optin} {
 
  upvar $optin options
  variable RTD
  foreach i { ID_REF_LrRNA ID_REF_SrRNA_1 ID_REF_SrRNA_2 ID_REF_SrRNA_3 ID_REF_SrRNA_4 ID_REF_SrRNA_5 COMMON_LrRNA COMMON_BODY COMMON_HEAD } { 
   set $i $RTD($i)
  }
  if { ! [info exists options(ONLYSSU)] } {
   drawtube $ID_REF_LrRNA $COMMON_LrRNA 1.8 0
  }
  if { ! [info exists options(ONLYLSU)] } {
   drawtube $ID_REF_SrRNA_1 $COMMON_BODY 1.8 1
   drawreptube $ID_REF_SrRNA_1 "all" 1.7 7
   drawtube $ID_REF_SrRNA_2 $COMMON_HEAD 1.8 9
   drawreptube $ID_REF_SrRNA_2 "all" 1.7 4
  }
  if { ! [info exists options(ONLYSSU)] && ! [info exists options(ONLYLSU)] } {
   drawtube $ID_REF_SrRNA_3 "all" 1.8 4
   mol off $ID_REF_SrRNA_3 
  }
  if { ! [info exists options(ONLYLSU)] } {
   drawtube $ID_REF_SrRNA_4 "all" 1.8 32
   mol off $ID_REF_SrRNA_4
 
   drawtube $ID_REF_SrRNA_5 "928 to 1389" 1.83 23
   drawreptube $ID_REF_SrRNA_5 "0 to 927 1390 to 1600" 1.83 35 "on"
   mol off $ID_REF_SrRNA_5
  }
 }
 
 proc findbothtrios {} {
 
  variable radenv
  variable RTD
  foreach i { ID_REF_SrRNA_1 ID_REF_SrRNA_2 ID_REF_SrRNA_3 ID_REF_SrRNA_4 COMMON_BODY COMMON_HEAD } { 
   set $i $RTD($i)
  }
  if { [info exists radenv(FINDHEADTRIO)] } {
   puts "Will try to find 3 residues that lie in a plane for the head"
   # the list at the end determines the + direction for the plane 
   findbesttrio [list $ID_REF_SrRNA_4 $ID_REF_SrRNA_2 ] $COMMON_HEAD 50 30 [list 1078 1136]
  }
  if { [info exists radenv(FINDBODYTRIO)] } {
   puts "Will try to find 3 residues that lie in a plane for the body"
   findbesttrio [list $ID_REF_SrRNA_3 $ID_REF_SrRNA_1] $COMMON_BODY 200 60 [list 697 843]
  }
 }
  
 proc defineresidues {} {
 
  # Define which residues will be used for vector calculations
  # These residues were chosen since they minimized the tilting value when measuring angles for 4V9D (classical and body rotated).  For the head, the atoms minimized the tilting value found when comparing 4V9D and 4V4Q
 
  # specific PDB chains  (large rRNA bundle, chain; small rRNA bundle, chain): 
 	# Body/head Unrotated 4V9D (3, E; 1, Y) 
 	# Body Rotated   4V9D (2, N; 1, A)
         # Head Rotated	 4V4Q (4, A; 3, A)
 
  # these were the values after 1/2022
 
  set body0 13 
  set body1 436
  set body2 823
  set body_d "21.638352083695878 83.94571620859452 25.344367696267405"
  set head0 1007
  set head1 1153 
  set head2 1315
  set head_d "65.11100908494743 8.627858902513841 57.50606947725621"
 
 # which residues are used to define the zero direction for the tilting axis
  # body zero-tilt direction is roughly parallel to h44
  set body_tilt0_0 "1426 1474"
  set body_tilt0_1 "1449 1454"
  # head zero-tilt direction is roughly parallel to mRNA
  set head_tilt0_0 693 
  set head_tilt0_1 519
 
  # to set an absolute origin for defining translation coordinates, we will use the P atoms of
  # resid $absorig in the 16S reference body model.  This is an arbitrarily chosen origin
  # and the absolute translations are only written with dump.
  set absorig 1409
  # To find the atom sets, the findbesttrio proc was used.  The last usage was on 1/1/22 (PCW).
  # Exact command used:
  #  BODY
  #   get the best trio (i.e. in a plane)
  #   set radenv(FINDBODYTRIO) 0
  #   RADtool -l 4v9d-pdb-bundle2.pdb -lc N -s 4v9d-pdb-bundle1.pdb -sc A  -stamp -findhead
  #   unset radenv(FINDBODYTRIO)
  #   Update body<N> values and then find the distances for triangulation
  #   set radenv(FIND_CENTER_ON) 0
  #   source RADtool.tcl
  #   RADtool -l 4v9d-pdb-bundle2.pdb -lc N -s 4v9d-pdb-bundle1.pdb -sc A  -stamp -findhead
  #   unset radenv(FIND_CENTER_ON)
  #   The distances returned are assigned to the variable body_d, below.
  # tilt value of 8.526471130583128e-5 was found
  
  #  HEAD 
  #   get the best trio (i.e. in a plane)
  #   set radenv(FINDHEADTRIO) 0
  #   RADtool -l 4v4q-pdb-bundle4.pdb -lc A -s 4v4q-pdb-bundle3.pdb -sc A  -stamp -findhead
  #   unset radenv(FINDHEADTRIO)
  #   Update head<N> values and then find the distances for triangulation
  #   set radenv(FIND_CENTER_ON) 0
  #   source RADtool.tcl
  #   RADtool -l 4v4q-pdb-bundle4.pdb -lc A -s 4v4q-pdb-bundle3.pdb -sc A  -stamp -findhead
  #   unset radenv(FIND_CENTER_ON)
  #   change values of head_d, below.
  # tilt value of 0.0005611904110883216 was found
  
  #####################OLD VALUES#########################
  # these were the values for 10/2021 to 1/2022
  #set body0 243  
  #set body1 810
  #set body2 443
  #set body_d "15.356563966567302 30.222420991707867 76.70248662566449"
  #set head0 1065 
  #set head1 1215
  #set head2 1336 
  #set head_d "18.53213084091104 57.219988731536354 57.020821095036524"
  
  # These were the resids and distances used before 10/2021
  # set body0 41 
  # set body1 127 
  # set body2 911
  # set body_d "64.27260243567858 52.429574225326476 19.250916881084546"
  #
  # set head0 984 
  # set head1 940 
  # set head2 1106
  # set head_d "12.420821202711673 41.506695819235986 55.19369134425787"
  
  #################END OLD VALUES#########################
  AddToRTD body0 body1 body2 head0 head1 head2 body_tilt0_0 body_tilt0_1 head_tilt0_0 head_tilt0_1 body_d head_d absorig 
 }
 
 proc loadmodel {optin} {
  variable fullfilenames 
  upvar $optin options
  set trajlist ""
  puts -nonewline $options(OUTPUT) "\nStructural model/s of the "
  if { [info exists options(ONLYSSU)] } {
   puts $options(OUTPUT) "small subunit"
  } elseif { [info exists options(ONLYLSU)] } {
   puts $options(OUTPUT) "large subunit"
  } else {
   puts $options(OUTPUT) "large and small subunits"
  }
 
  if { [info exists options(LOADTRAJ)] } {
   set trajlist $options(LOADTRAJ)
   puts $options(OUTPUT) "Will load a trajectory stored in the file/s $options(LOADTRAJ)"
   if {[catch {
    set nm [lindex $trajlist 0]
    set ID_LARGE  [mol new $nm waitfor -1]
    set fullfilenames($ID_LARGE) [file normalize $nm]
   } errorMessage] != 0} {
    error "Unable to load [lindex $trajlist 0]. Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
   }
   set ID_SMALL $ID_LARGE
  } else {
   if { ! [info exists options(ONLYSSU)] } {
    if { [info exists options(LARGEPDB)] } {
     set LARGEPDB $options(LARGEPDB)
     puts $options(OUTPUT) "    LSU rRNA structure file: $LARGEPDB"
     if {[catch {
      set ID_LARGE  [mol new $LARGEPDB waitfor -1]
      set fullfilenames($ID_LARGE) [file normalize $LARGEPDB]
     } errorMessage] != 0} {
      error "Unable to load $LARGEPDB. Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
     }
    } 
   } else {
    # since we are only looking at the head, we will set LSU stuff to null
    set ID_LARGE -1
    set CHAIN_ID_LrRNA -1
   } 
   if { [info exists options(SMALLPDB)] } {
    set SMALLPDB $options(SMALLPDB)
    puts $options(OUTPUT) "    SSU rRNA structure file: $SMALLPDB"
    if { [info exists options(LARGEPDB)] && $LARGEPDB eq $SMALLPDB } {
     set ID_SMALL $ID_LARGE
    } else {
     if {[catch {
      set ID_SMALL  [mol new $SMALLPDB waitfor -1] 
      set fullfilenames($ID_SMALL) [file normalize $SMALLPDB]
     } errorMessage] != 0} {
      error "Unable to load $SMALLPDB. Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
     }
    }
   } else {
    set ID_SMALL -1
    set CHAIN_ID_SrRNA -1
   } 
  }
 
  if { ! [info exists options(ONLYSSU)] } {
   if { [info exists options(CHAIN_LrRNA)] } {
    set CHAIN_ID_LrRNA $options(CHAIN_LrRNA)
    puts $options(OUTPUT) "    large rRNA chain ID: \"$CHAIN_ID_LrRNA\""
   } 
  }
 
  if { ! [info exists options(ONLYLSU)] } {
   if { [info exists options(CHAIN_SrRNA)] } {
    set CHAIN_ID_SrRNA $options(CHAIN_SrRNA) 
    puts $options(OUTPUT) "    small rRNA chain ID: \"$CHAIN_ID_SrRNA\""
   } 
 
   if { [info exists options(STAMP)] } {
    if {[info exists options(HEADFIRST)]} {
     puts $options(OUTPUT) "    STAMP alignment info: small rRNA, first head residue considered: $options(HEADFIRST)"
    }
    if {[info exists options(HEADLAST)]} {
     puts $options(OUTPUT) "    STAMP alignment info: small rRNA, last head residue considered: $options(HEADLAST)"
    }
   }
  }
  AddToRTD ID_LARGE CHAIN_ID_LrRNA ID_SMALL CHAIN_ID_SrRNA 
 }
 
 proc drawreptubelines {MOL_ID CHAIN_ID WIDTH COLOR {addnewrep 1}} {
 
  if { $addnewrep == 0 } {
   mol addrep $MOL_ID
  }
  # color the first rep, which should be lines.
  set repid [expr [molinfo $MOL_ID get numreps] - 1] 
  mol modselect $repid $MOL_ID chain \"$CHAIN_ID\" 
  mol modcolor $repid $MOL_ID ColorID $COLOR 
 
  # make a rep and draw in tubes
 
  mol addrep $MOL_ID
  set repid [expr [molinfo $MOL_ID get numreps] - 1] 
  mol modselect $repid $MOL_ID chain \"$CHAIN_ID\" 
  mol modstyle $repid $MOL_ID Tube $WIDTH 12.000000
  mol modcolor $repid $MOL_ID ColorID $COLOR 
 }
 
 proc drawreptube {MOL_ID RESIDUES WIDTH COLOR {onoff "off"}} {
 
  # make a rep and draw in tubes
 
  mol addrep $MOL_ID
  set repid [expr [molinfo $MOL_ID get numreps] - 1] 
  if { $RESIDUES == "all" } {
   mol modselect $repid $MOL_ID "all" 
  } elseif { [regexp -nocase {^(\s+)?[0-9]} $RESIDUES ] != 1  } {
   mol modselect $repid $MOL_ID "$RESIDUES" 
  } else {
   mol modselect $repid $MOL_ID "resid $RESIDUES" 
  }
  mol modstyle $repid $MOL_ID Tube $WIDTH 12.000000
  mol modcolor $repid $MOL_ID ColorID $COLOR 
  mol showrep $MOL_ID $repid $onoff
 }
 
 proc drawtube {ID RESIDUES WIDTH COLOR} {
 
  if { $RESIDUES == "all" } {
   mol modselect 0 $ID "all"
  } elseif { [regexp -nocase {^(\s+)?[0-9]} $RESIDUES ] != 1  } {
   mol modselect 0 $ID "$RESIDUES"
  } else {
   mol modselect 0 $ID "resid $RESIDUES"
  }
  mol modstyle  0 $ID Tube $WIDTH 12.000000
  mol modcolor  0 $ID ColorID $COLOR
 }
 
 proc calc_body_angle { optin p0_LrRNA p1_LrRNA p1_body p0_body zero_body large small framen ID_LARGE ID_SMALL ID_REF_SrRNA_1 ID_REF_SrRNA_4 ID_REF_SrRNA_5 ref_norv_body ref_perp_body ref_bond_body body0 body1 body2 body_d body_c_0 ref16_1 ref16_2 head_tilt0_0 head_tilt0_1 body_absorig } {
 
  variable radenv
  upvar $optin options
  # Calculates the body rotation angle
 
  if { ! [info exists options(ONLYSSU)] } {
   ####### Align both the large and small subunits to the reference large rRNA
   set M [measure fit $p1_LrRNA $p0_LrRNA]
  #Report Mlarge:
   set Mlarge $M
 
   $large frame $framen
   $large move $M
   set RMSD_LrRNA [measure rmsd $p0_LrRNA $p1_LrRNA]
   if { [info exists options(ONLYLSU)] } {
    return $RMSD_LrRNA
   }
   if { $ID_LARGE != $ID_SMALL } {
    $small frame $framen
    $small move $M
   }
  } else {
   set RMSD_LrRNA -1
  }
 
  if { [info exists options(isTraj)] } {
   set body_ref_ax_pos "not defined"
  } else {
   if { [info exists radenv(FIND_CENTER_ON)] } {
    # draw the vector for the reference body orientation
    set body_ref_ax_pos {1 1 1}
   } else {
    set body_ref_ax_pos [get_rot_center1 $ID_REF_SrRNA_1 $ref_norv_body $body0 $body1 $body2 $body_d ]
   }
  }
 
  ###### Idealize the coordinates of the SSU body
  if { ! [info exists options(ONLYSSU)] } {
   set M [measure fit $p0_body $p1_body]
  # Report Mbody:
   set Mbody $M
 
   $ref16_1 move $M
   $ref16_2 move $M
   set mysel [atomselect $ID_REF_SrRNA_4 "all"]
   $mysel move $M
   $mysel delete
  } else {
   # we are not only looking at the head, then don't move the 30S
   # if we are only looking at the head, then align to the reference body
   set M [measure fit $p1_body $p0_body]
  # Report Mbody:
   set Mbody $M
 
   $small move $M
  }
  set RMSD_body [measure rmsd $p1_body $p0_body]
  # TODO before fitting the head, we should get the zero direction of the head tilt
  set zero_head [ get_one_vector $ID_REF_SrRNA_1 $head_tilt0_0 $head_tilt0_1 ] 
 
  if { [info exists options(ONLYSSU)] } {
   return [list $zero_head $RMSD_body $Mbody {0 0 0}]
  }
  ###### Get the vectors "bond" and "norv" of the rotated rigid SrRNA 		
  lassign [get_vectors  $ID_REF_SrRNA_1 $body0 $body1 $body2]  bond norv perp
 
  if { [info exists radenv(FIND_CENTER_ON)] } {
   # if we are going to find the center of rotation, keep track of an atom
   set body_c_1 [get_xyz $ID_REF_SrRNA_1 $body1]  
   # this def of the ax pos is to avoid an error in get_rot_center1
   # get_rot_center1 uses triangulation to find the position of the axis. 
   # But, FIND_CENTER_ON means we are trying to find its position.  
   set body_fit_ax_pos {0 0 0}
  } else {
   if { ! [info exists options(isTraj)] } {
    set body_fit_ax_pos [get_rot_center1 $ID_REF_SrRNA_1 $norv $body0 $body1 $body2 $body_d ]
   } else {
    set body_fit_ax_pos {0 0 0}
   } 
  }
 
  ###### Calculate the rotations of the SSU body
  lassign [calc_euler $ref_norv_body $norv $zero_body $ref_bond_body $bond "body" ] body_rotation body_phi body_psi body_tilt body_tilt_direction body_tilt_vector zero_body bodytvecs  
  lassign [find_ER_between_2_systems $ref_bond_body $ref_perp_body $ref_norv_body $bond $perp $norv ] ERangle_body ERvec_body
  
  if { [info exists options(isTraj)] } {
   # not currently supported with trajectories
   set deltaabs_body {0 0 0}
  } else {
   set deltaabs_body [transdiff $body_absorig $body_ref_ax_pos $body_fit_ax_pos $ERangle_body $ERvec_body ]
  }
  if { ! [info exists options(isTraj)]  } {
 
   foreach name {body_ref_ax_pos body_fit_ax_pos norv ref_norv_body bond ref_bond_body perp ref_perp_body ERvec_body } {
    set j [set "$name"]
    if {$j ne "not defined"} {
     set $name [format_list $j "%.8f"]
    }
   }
 
   set axes_body_ref "X: $body_ref_ax_pos\n            V1: $ref_bond_body\n            V2: $ref_perp_body\n            V3: $ref_norv_body"
   set axes_body_fit "X: $body_fit_ax_pos\n            V1: $bond\n            V2: $perp\n            V3: $norv"
  } else {
 
   set axes_body_ref "not defined" 
   set axes_body_fit "not defined"
  }
 
  ## if we are finding the centers of rotation, then do the following...
  if { [info exists radenv(FIND_CENTER_ON)] } {
   set start [calc_rot_center $ref_norv_body $body_rotation $body_c_0 $body_c_1 $options(OUTPUT)] 
   set xyz0 [get_xyz $ID_REF_SrRNA_1 $body0]
   set xyz1 [get_xyz $ID_REF_SrRNA_1 $body1]
   set xyz2 [get_xyz $ID_REF_SrRNA_1 $body2]
   set dist0 [veclength [vecsub $start $xyz0]]
   set dist1 [veclength [vecsub $start $xyz1]]
   set dist2 [veclength [vecsub $start $xyz2]]
   puts $options(OUTPUT) "\n\nFINDING AXES OUTPUT: BODY"
   puts $options(OUTPUT) "\tPlane defined by ResiDs $body0 $body1 $body2"
   puts $options(OUTPUT) "\tPosition of axis defined by distances:$dist0 $dist1 $dist2\n"
  }
  return [list $zero_head $RMSD_LrRNA $RMSD_body $body_rotation $body_phi $body_psi $body_tilt $body_tilt_direction $body_tilt_vector $axes_body_ref $axes_body_fit $ERvec_body $ERangle_body $body_ref_ax_pos $ref_norv_body $body_fit_ax_pos $norv $zero_body $Mlarge $Mbody $deltaabs_body $bodytvecs]
 }
 
 proc calc_head_angle {optin p1_head p0_head zero_head ref16_2 ID_REF_SrRNA_1 ID_REF_SrRNA_2 head0 head1 head2 head_d head_absorig} {
 
  variable radenv
  upvar $optin options
 
  ####### Get the vectors "ref_bond" and "ref_norv" of the reference head
  lassign [get_vectors $ID_REF_SrRNA_2 $head0 $head1 $head2] ref_bond_head ref_norv_head ref_perp_head
  if { [info exists radenv(FIND_CENTER_ON)] } {
   # if we are going to find the center of rotation, keep track of an atom
   set head_c_0 [get_xyz $ID_REF_SrRNA_2 $head1]  
  }
 
  if { [info exists options(isTraj)]  } {
   set head_ref_ax_pos "not defined"
  } else {
  ###### Get the vectors "bond" and "norv" of the rotated rigid SrRNA 		
   if { [info exists radenv(FIND_CENTER_ON)] } {
    # draw the vector for the reference body orientation
    set head_ref_ax_pos {1 1 1}
   } else {
    set head_ref_ax_pos [get_rot_center1 $ID_REF_SrRNA_2 $ref_norv_head $head0 $head1 $head2 $head_d ]
   }
  }
  ###### Idealize the coordinates of the SSU head
  set M [measure fit $p0_head $p1_head]
  # we need both M and M2, so we can undo the alignment of p0_head, back to the body reference
  set M2 [measure fit $p1_head $p0_head]
  # Report Mhead:
  set Mhead $M
 
  $ref16_2 move $M
  set RMSD_head [measure rmsd $p1_head $p0_head]
  ###### Get the vectors "bond" and "norv" of the rotated rigid head
  lassign [get_vectors $ID_REF_SrRNA_2 $head0 $head1 $head2] bond norv perp 
 
  if { ! [info exists options(isTraj)]  } {
   # draw the vector for the rotated head 
   if { [info exists radenv(FIND_CENTER_ON)] } {
    # draw the vector for the reference body orientation
    set head_fit_ax_pos {0 0 0}
   } else {
    set head_fit_ax_pos [get_rot_center1 $ID_REF_SrRNA_2 $norv $head0 $head1 $head2 $head_d ]
   }
  } else {
   set head_fit_ax_pos {0 0 0}
  }
 
  if { [info exists radenv(FIND_CENTER_ON)] } {
   # if we are going to find the center of rotation, keep track of an atom
   set head_c_1 [get_xyz $ID_REF_SrRNA_2 $head1]
  }
 
  ###### Calculate the rotations of the SSU head
  lassign [calc_euler $ref_norv_head $norv $zero_head $ref_bond_head $bond "head" ] head_rotation head_phi head_psi head_tilt head_tilt_direction head_tilt_vector zero_head headtvecs 
  lassign [find_ER_between_2_systems $ref_bond_head $ref_perp_head $ref_norv_head $bond $perp $norv ] ERangle_head ERvec_head
  if { [info exists options(isTraj)] } {
   set deltaabs_head {0 0 0}
  } else {
   set deltaabs_head [transdiff $head_absorig $head_ref_ax_pos $head_fit_ax_pos $ERangle_head $ERvec_head ]
  }
  if { ! [info exists options(isTraj)]  } {
 
   foreach name {head_ref_ax_pos head_fit_ax_pos norv ref_norv_head bond ref_bond_head perp ref_perp_head ERvec_head } {
    set j [set $name]
    if {$j ne "not defined"} {
     set $name [format_list $j "%.8f"]
    }
   }
 
   set axes_head_ref "X: $head_ref_ax_pos\n            V1: $ref_bond_head\n            V2: $ref_perp_head\n            V3: $ref_norv_head"
   set axes_head_fit "X: $head_fit_ax_pos\n            V1: $bond\n            V2: $perp\n            V3: $norv"
 
  } else {
 
   set axes_head_ref "not defined" 
   set axes_head_fit "not defined"
  }
 
  ## if we are finding the centers of rotation, then do the following...
  if { [info exists radenv(FIND_CENTER_ON)] } {
   set start [calc_rot_center $ref_norv_head $head_rotation $head_c_0 $head_c_1 $options(OUTPUT)] 
 
  if { ! [info exists options(isTraj)]  } {
   set tmp [mol new]
   draw_arrow_scaled $tmp $start $ref_norv_head 100 100 1
   mol rename $tmp "center of rotation: head"
  }
   set xyz0 [get_xyz $ID_REF_SrRNA_2 $head0]
   set xyz1 [get_xyz $ID_REF_SrRNA_2 $head1]
   set xyz2 [get_xyz $ID_REF_SrRNA_2 $head2]
   set dist0 [veclength [vecsub $start $xyz0]]
   set dist1 [veclength [vecsub $start $xyz1]]
   set dist2 [veclength [vecsub $start $xyz2]]
   puts $options(OUTPUT) "\n\nFINDING AXES OUTPUT: HEAD"
   puts $options(OUTPUT) "\tPlane defined by ResiDs $head0 $head1 $head2"
   puts $options(OUTPUT) "\tPosition of axis defined by distances:$dist0 $dist1 $dist2\n"
  }
  return [list $RMSD_head $head_rotation $head_phi $head_psi $head_tilt $head_tilt_direction $head_tilt_vector $axes_head_ref $axes_head_fit $ERvec_head $ERangle_head $head_ref_ax_pos $ref_norv_head $head_fit_ax_pos $norv $zero_head $Mhead $M2 $deltaabs_head $headtvecs]
 }
 
 proc dostampONLY {optin ID_TMP CHAIN_ID maph REFPDB refchain} {
 
  package require seqdata 1.1
  # ID_TMP - Mol ID of the model we are analyzing.
  # CHAIN_ID - chain ID of the model we are analyzing
  # maph - mapping file between RESIDs and residue numbers in model
  # REFPDB - option describing which reference model to use: LARGEREF, HEADREF or BODYREF
  # refchain - chain ID of the reference model. 
  # The main difference with this version of dostampONLY is that this version take a mapping file of the residue numbers, and also returns a list of residues (residue numbers) that are aligned.  This is useful for multiple iterations of STAMP. Also, this routine only works on a single regions at a time, not all three regions.
  upvar $maph map
  upvar $optin options 
  set SCRIPTPATH $options(SCRIPTPATH)
 
  set ID_ALIGN  [mol new "$SCRIPTPATH/share/reference_models/$REFPDB" waitfor -1] 
 
  # set the selections for the reference E. coli sequence
  set selt [atomselect  $ID_ALIGN "name P"]
  set sel_reslist [$selt get residue]
  $selt delete
  # set the selections for the model that we are analyzing. Note, we already pre-processed the selections, so now we just select "all"
  set selt [atomselect $ID_TMP "all"]
  set sel_model_reslist [lsort -unique -integer [$selt get residue]]
  $selt delete
 
  # align and save the sequence ID of the aligned sequences
  lassign [::RTALIGN::single_stamp $options(TMPNAME) $ID_TMP $CHAIN_ID $sel_model_reslist $ID_ALIGN "$refchain" $sel_reslist] seqIDmod seqIDalign
  # store the corresponding resids, so that we can change the model IDs before proceeding with angle calculations
  set nalign 0
  set newlist ""
  #make a mapping array from original numbers to E coli-aligned numbering
  for {set i 0} { $i <[::SeqData::getSeqLength $seqIDalign]} {incr i} {
   set r1 [::SeqData::getResidueForElement $seqIDalign $i]
   set r2 [::SeqData::getResidueForElement $seqIDmod $i]
   set n1 [::SeqData::getElement $seqIDalign $i]
   set n2 [::SeqData::getElement $seqIDmod $i]
   lassign [cleanup $r1] r1 ins1
   lassign [cleanup $r2] r2 ins2
 
   if { $n1 ne "-" && $n2 ne "-"} {
    incr nalign
    if {![info exists newmap($r2$ins2)]} {
     set newmap($r2$ins2) $r1$ins1
     lappend newlist $map($r2$ins2)
    } else {
     error "Duplicated aligned residue ($r2$ins2) in small rRNA. step 2"
    }
   }
  }
 
  mol delete $ID_ALIGN
  return [list $nalign [array get newmap] $newlist]
 
 }
 
 proc doseqstampalign { optin } {
  variable radversion
  variable fullfilenames 
  upvar $optin options
  set headstart 928 
  set headend 1389 
  if { [info exists options(FINDHEAD)] } {
   set seqlistonly 1
  } else {
   set seqlistonly 0
   # in findseqaligned, only get a list of residue numbers, don't perform sequence alignment.
  }
  
  set tind [findseqaligned options $seqlistonly] 
 
  if { ! [info exists options(ONLYLSU)] } {
   set CHAIN_ID_SrRNA $options(CHAIN_SrRNA)
   array set smallmap [lindex $tind 1]
   # by default, don't split the small subunit for alignment
   set sephead "all"
   set sepbody "all"
   set h_head_b 0
   set b_head_b 0
   set h_head_e 0
   set b_head_e 0
   if {[info exists options(HEADFIRST)] || [info exists options(HEADLAST)]} {
    # since we gave definitions for the head, separately align it during STAMP
    set sephead "head"
    set sepbody "body"
    set h_head_b $options(HEADFIRST)
    set h_head_e $options(HEADLAST)
    set b_head_b $options(HEADFIRST)
    set b_head_e $options(HEADLAST)
   }
   if {$seqlistonly == 1} {
    # guess where the head is...
    # 
    array set smap [lindex $tind 2 ]
    set nearF 1000000
    set nearL 1000000
    foreach name [array names smap] {
     if {[expr $smap($name)-900] > 0 } {
      set diff [expr $smap($name)-900]
     } else {
      set diff [expr 900-$smap($name)]
     } 
     if {$diff < $nearF } {
      set nearF $diff
      set closestB $name
     }
    }
    if { ! [info exists closestB] } {
     set arsize [array size smap]
     error "Unable to find a suitable guess for where the SSU head may start. This often means sequence alignment failed. Sequence alignment aligned $arsize residues."
    }
    if {[catch {
     set ID_SMALL [mol new "$options(SMALLPDB)" waitfor -1]
     set fullfilenames($ID_SMALL) [file normalize $options(SMALLPDB)]
    } errorMessage] != 0} {
     error "Unable to load $options(SMALLPDB). Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
    }
    # convert to residue numbering, in case there are issues with the resids (e.g. accidental jumps in sequence)
    set first [atomselect $ID_SMALL "name P and chain \"$CHAIN_ID_SrRNA\" and resid $closestB"]
    set firstresidue [$first get residue]
    $first delete
   
    if {[info exists options(FINDHEAD2)]} { 
     set resrange1 "residue $firstresidue to [expr $firstresidue+200]"
     set resrange2 "residue [expr $firstresidue+200] to [expr $firstresidue+500]"
    } else {
     set resrange1 "resid $closestB to [expr $closestB+200]"
     set resrange2 "resid [expr $closestB+200] to [expr $closestB+500]"
    }

    set tentativehead1 [atomselect $ID_SMALL "not backbone and not name P and chain \"$CHAIN_ID_SrRNA\" and $resrange1"]
    set tentativehead2 [atomselect $ID_SMALL "not backbone and not name P and chain \"$CHAIN_ID_SrRNA\" and $resrange2"]
    set dthresh 3
    if { [$tentativehead1 num] == 0 || [$tentativehead2 num] == 0 } {
     # this usually only happens if we have a P-only model.  so, we will just look at P atoms, instead, and use a longer threshhold
     set tentativehead1 [atomselect $ID_SMALL "name P and chain \"$CHAIN_ID_SrRNA\" and $resrange1"]
     set tentativehead2 [atomselect $ID_SMALL "name P and chain \"$CHAIN_ID_SrRNA\" and $resrange2"]
     set dthresh 24
    }

    set conts [ measure contacts $dthresh $tentativehead1 $tentativehead2]
    set nconts [llength [lindex $conts 0]]
    if { $nconts == 0 } {
     error "sequence-based alignment could not make a good guess about where the SSU head is located. This error is common when trying to analyze atypical ribosomes (e.g. those that are composed of many rRNA molecules per subunit). If this is a structure composed of large rRNA molecules in the LSU and SSU (e.g. 23S and 16S), then perhaps try using the -h_f and -h_l options. It can also happen if there are issues in the head residue numbering. The flag -findhead2 (Alt. find head method button) can be helpful." 
    }
    $tentativehead1 delete
    $tentativehead2 delete
    set diff 0
    set c1 [lindex $conts 0]
    set c2 [lindex $conts 1]
    for {set i 0} { $i < $nconts} {incr i} {
     set a1 [lindex $c1 $i] 
     set a2 [lindex $c2 $i]
     set dt [expr $a2-$a1]
     if { $dt > $diff } {
      set at1 $a1
      set at2 $a2
      set diff $dt
     }
    }
    set selt [atomselect $ID_SMALL "index $at1"]
    set closestB [$selt get resid]
    $selt delete
    set selt [atomselect $ID_SMALL "index $at2"]
    set closestE [$selt get resid]
    $selt delete
    mol delete $ID_SMALL
    set headstart $closestB 
    set headend $closestE
 
    puts $options(OUTPUT) "When performing STAMP alignment, will try to align E. coli head to residues $closestB to $closestE of the model"
 
    set h_head_b $closestB
    set h_head_e $closestE
    set b_head_b $closestB
    set b_head_e $closestE
    set options(HEADFIRST) $closestB
    set options(HEADLAST)  $closestE
    set sephead "head"
    set sepbody "body"
   }
 
   puts $options(OUTPUT) "\nstarting STAMP alignment of small subunit rRNA, head"
   array set stamphead  [dostructurealignments options $options(SMALLPDB) $CHAIN_ID_SrRNA smallmap $options(HEADREF) "Y" $sephead $h_head_b $h_head_e]
   puts $options(OUTPUT) "\nstarting STAMP alignment of small subunit rRNA, body"
   array set stampbody [dostructurealignments options  $options(SMALLPDB) $CHAIN_ID_SrRNA smallmap $options(BODYREF) "Y" $sepbody $b_head_b $b_head_e]
  }
 
  if { ! [info exists options(ONLYSSU)] } {
   set CHAIN_ID_LrRNA $options(CHAIN_LrRNA)
   array set largemap [lindex $tind 0]
   puts $options(OUTPUT) "\nstarting STAMP alignment of large subunit rRNA"
   array set stamplarge [dostructurealignments options $options(LARGEPDB) $CHAIN_ID_LrRNA largemap $options(LARGEREF) "E" "all" 0 0]
  } else {
   set stamplarge(0) 0
  }
 
  if { ! [info exists options(ONLYLSU)] } {
   foreach name [array names stamphead] {
    set finalsmall($name) $stamphead($name)
   }
   foreach name [array names stampbody] {
    if {[info exists finalsmall($name)]} {
     if { $finalsmall($name) ne $stampbody($name) } {
      # hope we never hit this error
      error "During STAMP alignment, a residue was found to align in the head and body, but the aligned numbers are not the same.  Issue associated with residues $name and $finalsmall($name)"
     }
    }
    set finalsmall($name) $stampbody($name)
   }
  } else {
   set finalsmall(0) 0
  }
 
  if {[info exists options(ALIGNOUT)] } {
   puts $options(ALIGNOUT) "Mapping to Ecoli sequencing used for angle calculations.  The resid in the input model is followed by the E. coli number.
Generated by RADtool v$radversion, on [clock format [clock second] -format "%D %T"], using machine [info hostname]"
   # write out final mapping:
  
   if { ! [info exists options(ONLYSSU)] } {
    puts -nonewline $options(ALIGNOUT) "large subunit: "
    foreach name [lsort  -dictionary  [array names stamplarge]] {
     puts -nonewline $options(ALIGNOUT)  "$name $stamplarge($name) "
    }
    puts -nonewline $options(ALIGNOUT)  "\n"
   }
 
   if { ! [info exists options(ONLYLSU)] } {
    puts -nonewline $options(ALIGNOUT) "small subunit: "
    foreach name [lsort -dictionary [array names finalsmall]] {
     puts -nonewline $options(ALIGNOUT) "$name $finalsmall($name) "
    }
    puts -nonewline $options(ALIGNOUT)  "\n"
    puts  $options(ALIGNOUT) "small head: $headstart to $headend"
   }
   close $options(ALIGNOUT)
  }
  return [list [array get options] [array get stamplarge] [array get finalsmall]]
 }
 
 proc dostructurealignments {optin PDBNAME CHAIN_ID map PDBNAME_REF CHAIN_ID_REF lhb head_b head_e} {
 
  variable radenv
  variable fullfilenames 
  # PDBNAME - name of model.  e.g. options(SMALLPDB)
  # CHAIN_ID - model chain if
  # map - e.g. smallmap (resid to residue)
  # PDBNAME_REF e.g. options(HEADREF)
  # CHAIN_ID_REF .e.g. Y
  # also select based on the ranges provided.
  # save the file
  upvar $optin options
  upvar $map modmap
  set nalign 0
  set nalignlast -1
  # current set of aligned residues
  set lastlist ""
  set iter 0
  set lastnalign -1
  set nalign 0
  if {[catch {
   set ID_TMP      [mol new "$PDBNAME" waitfor -1] 
  } errorMessage] != 0} {
   error "Unable to load $PDBNAME. Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
  }
  # make a selection that includes the given chain ID, residues of the head that is also aligned before
  if { $lhb eq "head" } {
   set selstring "chain \"$CHAIN_ID\" and resid $head_b to $head_e" 
  } elseif { $lhb eq "body" } {
   set selstring "chain \"$CHAIN_ID\"  and (not resid $head_b to $head_e)"
  } elseif { $lhb eq "all" } { 
   set selstring "chain \"$CHAIN_ID\" " 
  } else {
   error "dostructurealignments error: invalid domain ref"
  }
 
  # we need to ignore alternate locators, since they break STAMP with nonsense error messages.
  # If there are only a few altlocs, this should not be a problem.
  set selstring "$selstring and name P and $radenv(NOTSTUFF) and not altloc B to Z"
 
  set tmpsel [atomselect $ID_TMP "$selstring"]
  if { [$tmpsel num] == 0 } {
   error "dostructurealignments error: 0 atoms identified by selection \"$selstring\""
  }
  $tmpsel set altloc {" "}
 
  if { [regexp -all "\.pdb$" $PDBNAME] == 1 } {
   set sname "$options(TMPNAME).pdb"
   animate write pdb "$sname" beg 0 end 0 skip 1 sel $tmpsel 
  } elseif { [regexp -all "\.cif$" $PDBNAME] == 1 } {
   set sname "$options(TMPNAME).cif"
   animate write pdbx "$sname" beg 0 end 0 skip 1 sel $tmpsel 
  }
 
  $tmpsel delete
  mol delete $ID_TMP
  set ID_TMP  [mol new "$sname" waitfor -1]
 
  set indt [dostampONLY options $ID_TMP $CHAIN_ID modmap $PDBNAME_REF "$CHAIN_ID_REF"]
  mol delete $ID_TMP
  catch {file delete $sname} tv 
  # the numbers that aligned after stamp
  set nalignc [lindex $indt 0]
  puts $options(OUTPUT) "$nalignc aligned residues"
 
  # the newmapping, based on stamp
  array set newmap [lindex $indt 1]
  # clear stored sequences
  ::SeqData::reset
  return [array get newmap]
 }
 
 proc findseqaligned { optin seqlistonly} {
  variable fullfilenames 
  upvar $optin options
 
  set SCRIPTPATH $options(SCRIPTPATH)
 
  if { ! [info exists options(ONLYSSU)] } {
   set CHAIN_ID_LrRNA $options(CHAIN_LrRNA)
   if {[catch {
    set ID_LARGE [mol new "$options(LARGEPDB)" waitfor -1]
    set fullfilenames($ID_LARGE) [file normalize $options(LARGEPDB)]
   } errorMessage] != 0} {
    error "Unable to load $options(LARGEPDB). Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
   }
   set selt [atomselect $ID_LARGE "chain \"$CHAIN_ID_LrRNA\" and name P and not altloc B to Z"]
   foreach nindex [$selt get index]  {
    set restmp [atomselect $ID_LARGE "index $nindex" ]
    set ins [$restmp get insertion]
    set residue [$restmp get residue]
    set resid [$restmp get resid]
    lassign [cleanup "$resid $ins"] resid ins
    set largemap($resid$ins) $residue
    $restmp delete
   }
   $selt delete 
   mol delete $ID_LARGE
 
  } else {
   # not used, so just give a summary entry
   set largemap(0) 0
  }
 
  if { ! [info exists options(ONLYLSU)] } {
   set CHAIN_ID_SrRNA $options(CHAIN_SrRNA)
 
   if {[info exists options(LOADTRAJ)]} {
    set SMALLPDB [ lindex $options(LOADTRAJ) 0]
   } else {
    set SMALLPDB $options(SMALLPDB)
   }
 
   if {[catch {
    set ID_SMALL [mol new "$SMALLPDB" waitfor -1]
    set fullfilenames($ID_SMALL) [file normalize $SMALLPDB]
   } errorMessage] != 0} {
    error "Unable to load $SMALLPDB. Something about the file is incompatible with VMD. Perhaps consult VMD documentation.\n\nError message:\n$errorMessage"
   }
   set selt [atomselect $ID_SMALL "chain \"$CHAIN_ID_SrRNA\" and name P and not altloc B to Z"]
   foreach nindex [$selt get index]  {
    set restmp [atomselect $ID_SMALL "index $nindex" ]
    set ins [$restmp get insertion]
    set residue [$restmp get residue]
    set resid [$restmp get resid]
    lassign [cleanup "$resid $ins"] resid ins
    set smallmap($resid$ins) $residue
    $restmp delete
   }
   $selt delete
  } else {
   # not used, so just give a summary entry
   set smallmap(0) 0
  }
  
  if { $seqlistonly == 0 } {
   # do not do sequence alignment.  Just save the mapping to residues
   mol delete $ID_SMALL
   return [list [array get largemap] [array get smallmap]]
  }
 
  set largemap(0) 0
 
  if { ! [info exists options(ONLYLSU)] } {
   # do SSU
   set ID_ALIGN_SrRNA [mol new "$SCRIPTPATH/share/reference_models/$options(SMALLREF)" waitfor -1] 
   set selt [atomselect  $ID_ALIGN_SrRNA "name P"]
   set sel_SrRNA_reslist [$selt get residue]
   $selt delete
   set selt [atomselect $ID_SMALL "chain \"$CHAIN_ID_SrRNA\" and name P and not altloc B to Z"]
   set sel_SrRNAmodel_reslist [$selt get residue]
   $selt delete
   if { [llength $sel_SrRNAmodel_reslist] == 0 } {
    error "Sequence Alignment Error: No atoms in small subunit selection \" chain $CHAIN_ID_SrRNA and name P\" "
   } 
   puts $options(OUTPUT) "\nstarting ClustalW sequence alignment of small subunit rRNA"
   lassign [::RTALIGN::single_seqalign $options(TMPNAME) $ID_SMALL $CHAIN_ID_SrRNA $sel_SrRNAmodel_reslist $ID_ALIGN_SrRNA "Y" $sel_SrRNA_reslist ] seqIDSrRNAmod seqIDSrRNAalign
   set alignlength [::SeqData::getSeqLength $seqIDSrRNAalign]
   if { $alignlength == 0 } {
    error "Alignment failed.  No residues aligned."
   } 
   set nalignSrRNA 0
   #make a mapping array from original numbers to E coli-aligned numbering
   for {set i 0} { $i < $alignlength} {incr i} {
    set r1 [::SeqData::getResidueForElement $seqIDSrRNAalign $i]
    set r2 [::SeqData::getResidueForElement $seqIDSrRNAmod $i]
    set n1 [::SeqData::getElement $seqIDSrRNAalign $i]
    set n2 [::SeqData::getElement $seqIDSrRNAmod $i]
    lassign [cleanup $r1] r1 ins1
    lassign [cleanup $r2] r2 ins2
    if { $n1 ne "-" && $n2 ne "-"} {
     incr nalignSrRNA
     # if the residue is aligned.
     if {![info exists SrRNAmap($r2$ins2)]} {
      if {[info exists smallmap($r2$ins2)]} {
       set SrRNAmap($r2$ins2) $r1$ins1
      }
     } else {
      error "Duplicated aligned residue ($r2$ins2) in step 1"
     }
    }
   }
   mol delete $ID_SMALL
   mol delete $ID_ALIGN_SrRNA
  } else {
   set SrRNAmap(0) 0
  }
  return [list [array get largemap] [array get smallmap] [array get SrRNAmap] ]
 }
 
 proc prunecore {optin MOL_ID CHAIN_ID COMMON ID_REF} {
  # this is my attempt at a rapid-prune process.  Rather than remove one residue per iteration, remove anything over 1A and add anything under.  Then, realign based on those under and add/remove.  Continue until convergence is reached.
 
  ####### Get common residues in a domain
  set SEL [atomselect $MOL_ID "chain \"$CHAIN_ID\" and resid $COMMON and name P and not insertion A to Z and not altloc B to Z"]
  upvar $optin options
  set TMPNAME [expr {int(rand()*1000000)}]
  set tmpname "$TMPNAME.pdb"
  $SEL writepdb "$tmpname"
 
  set newID [mol new $tmpname]
  set fullmod [atomselect $newID "resid $COMMON"]
  set fullref [atomselect $ID_REF "name P and resid $COMMON"]
 
  set diff 1
  set threshold $options(PRUNEBY)
  set threshold [ expr $threshold*$threshold ] 
  set newlist $COMMON
  set round 0 
  while { $diff > 0 } { 
   incr round
   set oldlist $newlist
   unset newlist
   set newlist {} 
   set selref [atomselect $ID_REF "name P and resid $oldlist"]
   set selmod [atomselect $newID "resid $oldlist"]
   # align the copy of the model to the reference structure
   set M [measure fit $selmod $selref]
   $fullmod move $M
   # check the displacement of each atom
   foreach c0 [$fullref get {x y z}] c1 [$fullmod get {x y z}] indx [$fullmod get resid] {
    set dx [expr [lindex $c0 0] - [lindex $c1 0]]
    set dy [expr [lindex $c0 1] - [lindex $c1 1]]
    set dz [expr [lindex $c0 2] - [lindex $c1 2]]
    set displacement2 [expr $dx*$dx + $dy*$dy + $dz*$dz] 
    if { $displacement2 < $threshold } {
     lappend newlist $indx
    }
   }
   $selref delete
   $selmod delete
   if {[llength $newlist] == 0} {
    error "During pruning, all residues were discarded. This generally happens for \"unusual\" ribosomes. However, sometimes this just means you need to assign E coli numbering via the -stamp flag"
   }
   if { $oldlist == $newlist } {
    set diff 0
   } else {
    set diff 1
   }
  }
  $fullmod delete
  $fullref delete
  set selt [atomselect $newID "resid $newlist"]
  set nn [$selt num]
  $selt delete
  puts $options(OUTPUT) "    pruning converged on $nn residues after $round iterations"
  mol delete $newID 
  catch {file delete $TMPNAME.pdb} tv
  return [list $nn $newlist]
 }
 
 proc getcommon {NAME MOL_ID CHAIN_ID CORE MIN outstr} {
  ####### Get common residues in a domain
  set SEL [atomselect $MOL_ID "chain \"$CHAIN_ID\" and resid $CORE and name P and not insertion A to Z and not altloc B to Z"]
  set NUM_COMMON [$SEL num]
  set COMMON [$SEL get resid]
  if { $NUM_COMMON == 0 } {
   error "Unable to match any residues in the $NAME reference and structure file. Perhaps you gave the wrong chain ID."
  }
  if {$NUM_COMMON < $MIN} { 
   puts $outstr "\n\n\nWARNING: only able to match [$SEL num] residues in the $NAME reference and structure file. \n\n\n"
  }
  return [list $NUM_COMMON $COMMON]
 }
 
 proc getallcommon {NAME REFNAME MOL_ID CHAIN_ID RESRANGE MIN outstr} {
  ####### Get common residues in a domain
  set selt [atomselect $REFNAME "name P and resid $RESRANGE and not insertion A to Z and not altloc B to Z"]
  set listall [$selt get resid]
  $selt delete
  set SEL [atomselect $MOL_ID "chain \"$CHAIN_ID\" and resid $listall and name P and not insertion A to Z and not altloc B to Z"]
  set NUM_COMMON [$SEL num]
  set COMMON [$SEL get resid]
  if { $NUM_COMMON == 0 } {
   error "Unable to match any residues in the $NAME reference and structure file. Perhaps you gave the wrong chain ID."
  }
  if {$NUM_COMMON < $MIN} { 
   puts $outstr "\n\n\nWARNING: only able to match [$SEL num] residues in the $NAME reference and structure file. \n\n\n"
  }
  return [list $NUM_COMMON $COMMON]
 }
 
proc printusage {} {

 variable radversion
 variable NGUYENREF
 variable RADREF
 puts "

       More usage information:
       This tool describes subunit orientations in terms of Euler
       Angles, as first applied to describe ribosome simulations in:

$NGUYENREF

       Since then, there have been many enhancements added, which are
       described in:

$RADREF

       Features include:
       	- support for a single structure file containing both subunits
        - support for mmCIF files
       	- support for trajectories
       	- support for non-E.coli residue numbering
        - support for isolated SSUs
        - STAMP alignment options
        - automated download of ribosome structures
        - automated detection of LSU-SSU rRNA pairs, isolated SSUs
              and isolated LSUs
       	- graphical representations of fitted configs, axes and planes
       	- Euler-Rodrigues angles are calculated (4V9D defined as 0)
        - calculation of translational displacements
        - various animation options
       	- reference values changed: The atoms used for determining 
       	      the rotation planes have been updated, relative to 
       	      the above reference. These new atoms were selected 
       	      since they minimized the tilt angle identified when 
       	      comparing PDB entries: 
                  4V9D:DA,BA (unrotated) and 4V9D:CA,AA (rotated body)
       		  4V9D:DA,BA (unrotated) and 4V4Q:DB,CA (rotated head)
       	- rigid-body alignment now uses 4V9D:DA,BA

While the GUI is suitable for most purposes, one may also use the command-based interface by issuing the command \"RADtool\" in the TkConsole of VMD.  

When using the command-line interface, the following flags are supported:

   Structural Model Information:
       -l <file>         : file that contains the large subunit
                             pdb or cif format are supported
       -lc <int/char>    : chain ID of large rRNA (23S in bacteria)
                             one letter for PDB, or 1+ letters for cif 
       -s <file>         : file that contains the small subunit (-l and -s may be the same) 
       -sc <int/char>    : chain ID of small rRNA (16S in bacteria)
       -bundle <list>    : list of names of structure files (PDB/cif) that should be analyzed.  
                             This option will indicate that all files should be read. 
                             RADtool will attempt to identify all LSU-SSU pairs (or just SSUs, 
                             if -SSUonly is used, or just LSUs if -LSUonly is given), and then 
                             analyze each identified ribosome
       -download <name>  : name is a 4-letter RCSB accession code.  RADtool will attempt to 
                             download the bundle, individual PDB, or cif file.
                             -download automatically activated -bundle capabilities 

   Output Options:
       -o <file>         : angle output file (stdout, if not given)
       -e <file>         : error file (stderr, if not given)
       -ot <file>        : trajectory angle output file
       -align_out <file> : output file to store the alignments
       -cores_out <file> : output file to store the utilized core residues
       -savepdb <file>   : name of files to save aligned LSU and SSU rRNA
                             file will be named <file>_SSU.pdb and <file>_LSU.pdb
       -overwrite        : overwrite any existing files
       -precision \[1\]    : number of decimals to include for calculated angles
       -dump             : save lots of stuff to the output file

   Alignment options:
       -stamp            : (recommended) perform STAMP alignment of large and small 
                             rRNA (body and head) to determine E.coli numbering
       -align_in <file>  : read alignment to E.coli numbering from file
       -cores_in <file>  : read the definitions of the cores from file
                             This implies -notall
       -SSUonly          : only calculate angles for SSU head, relative to SSU body
       -LSUonly          : only align the LSU to the reference and quit
       -findhead2        : use alternate method for guessing the head location
                             This can be useful if the head can not be 
                             automatically identified

   Additional STAMP options:
       -h_f              : first residue of small subunit head (ResID in the input model)
       -h_l              : last residue of small subunit head (ResID in the input model)
                             note: -h_f and -h_l should only be necessary if alignment
                                   failed when using -stamp alone
       -notall           : only use predefined \"core\" residues, and do not try to use 
                             all residues (in LSU, SSU head and SSU body) during structural 
                             alignment steps.  
       -noprune          : do not prune the core groups 
       -pruneby \[2\]      : calculate orientation based only on residues (P atoms) 
                             that are within N Angstroms (default 2) in the reference 
                             E. coli model and the input model 

   Animation option:
       -animate <seq>    : depict the orientation in terms of rotation (r), tilt (t) 
                             and translation (tr). 
                             seq is the sequence of motions to show.
                             Supported values: r_t_tr, t_r_tr, tr_r_t, tr_t_r
       -animate2         : depict the orientation change between two ribosomes analyzed
                             Ribosome IDs of start/end are given interactively
                             Number of frames is given interactively
                             Interactively indicate which structure model is visualized 

   Trajectory-specific options:
       -traj <list>      : list of files containing the trajectory (e.g. \"conf.pdb traj.xtc\") 
       -t_first \[0\] \ \    : first frame to analyze
       -t_last	\[-1\] \ \    : last frame to analyze. -1 indicates all.
       -t_step \[1\]  \ \    : analyze every every Nth frame
       -t_seglength <int>: load and analyze N frames as a time

   Other options:
       -help             : see the list of options
       -moreinfo         : see a detailed description of RADtool,
                             including visualization and analysis features
       -test             : run a series of internal checks to verify RADtool
                             is working properly
       -testall          : run a series of internal checks to verify RADtool
                             is working properly, including trajectory tests
       -testtraj         : run a series of internal checks to verify RADtool
                             is working properly, only test trajectory options


       For more information about this tool, see:
            http://radtool.org

       Please send feedback to Paul C. Whitford (p.whitford@northeastern.edu)."

}
 
 proc get_one_vector {ID_SrRNA res0 res1} {
 
  ####### Get the unit vector pointing from res0 to res1 (P atoms)
  set SEL [atomselect $ID_SrRNA "resid $res0 and not insertion A to Z"]
  if {[$SEL num] == 0} {
   error "Unable to find resid $res0 in molecule $ID_SrRNA"
  }
  set P_atom0 [geom_center $SEL]
  set SEL [atomselect $ID_SrRNA "resid $res1 and not insertion A to Z"]
  if {[$SEL num] == 0} {
   $SEL delete
   error "unable to find resid $res1 in molecule $ID_SrRNA"
   return
  }
  set P_atom1 [geom_center $SEL]
  $SEL delete
  set bond [vecnorm [vecsub $P_atom1 $P_atom0]]
  return $bond
 }
 
 proc get_xyz {ID_SrRNA res0} {
 
  set SEL [atomselect $ID_SrRNA "name P and resid $res0 and not insertion A to Z"]
  if {[$SEL num] != 1} {
   error "Unable to find P atom of resid $res0 in molecule $ID_SrRNA"
  }
  set P_atom0 [join [$SEL get {x y z}]]
  $SEL delete
  return $P_atom0
 }
 
 proc get_vectors {ID_SrRNA res0 res1 res2} {
 
  ####### Get the vectors "bond" and "norv"
  set SEL [atomselect $ID_SrRNA "name P and resid $res0 and not insertion A to Z"]
  if {[$SEL num] != 1} {
   error "Unable to find P atom of resid $res0 in molecule $ID_SrRNA"
  }
  set P_atom0 [join [$SEL get {x y z}]]
  set SEL [atomselect $ID_SrRNA "name P and resid $res1 and not insertion A to Z"]
  if {[$SEL num] != 1} {
   error "unable to find P atom of resid $res1 in molecule $ID_SrRNA"
   #return
  }
  set P_atom1 [join [$SEL get {x y z}]]
  set SEL [atomselect $ID_SrRNA "name P and resid $res2 and not insertion A to Z"]
  if {[$SEL num] != 1} {
   error "unable to find P atom of resid $res2 in molecule $ID_SrRNA"
   #return
  }
  set P_atom2 [join [$SEL get {x y z}]]
  set bond [vecnorm [vecsub $P_atom2 $P_atom0]]
  set secv [vecnorm [vecsub $P_atom1 $P_atom0]]
  set norv [vecnorm [veccross $bond $secv]]
  set perp [vecnorm [veccross $norv $bond]]
 
  return [list $bond $norv $perp]
 }
 
 proc calc_rot_center {ref_norv angle coord0 coord1 outstr} {
 
  set angle [expr $angle*$pi/180]
  # This is a proc that is used when defining new coordinates
  # The routine assumes the normal vector is known for a ``pure'' rotation/tilt
  # angle is the rotation, already calculated, about the normal
  # coord0 is a coordinate in the defined rotation plane (unrotated configuration)
  # coord1 is the same relative coordinate in the rotated configuration
  # note that we are only looking for a rotation axis, and this can not
  # describe any translational changes between the coordinates.  So, we must
  # add the translation to coord0
  # the output is a point along the axis of rotation 
  set proj [vecdot [vecsub $coord1 $coord0] $ref_norv]
  set coord0 [vecadd $coord0 [vecscale $proj $ref_norv]]
  set dR2 [vecscale 0.5 [vecsub $coord1 $coord0]]
  set dR2l [veclength $dR2]
  set halfdR [vecadd $coord0 $dR2] 
  set disp [expr $dR2l*tan(($pi - $angle)/2)]
  set to_center [vecnorm [veccross $ref_norv $dR2]]
  set to_center [vecscale $disp $to_center]
  set center [ vecadd $halfdR $to_center]
  return $center 
 }
 
 proc draw_fitted_axis {center vec color name {length 100} {scale 100} } {
 
  set object [mol new]
  mol rename $object "$name"
  draw_arrow_scaled $object $center $vec $length $scale $color
  mol off $object
 }
 
 proc draw_plane {center norm xdir color name {wx 100} {wy 100}} {
 
  # center is the center of the rectangle
  # norm is the normal vector
  # xdir is a direction for defining "x"
  # name is a name
  # xw and wy are the widths.  x is the direction of zero
  # if zero is not perp to norm, then the projection is substracted
  set norm [vecnorm $norm]
  set xdir [vecnorm $xdir]
  set proj [vecdot $norm $xdir] 
  set xdir [vecnorm [vecsub $xdir [vecscale $proj $norm]]]
  set ydir [vecnorm [veccross $xdir $norm]]
  set corner1 [vecadd $center [vecscale $wx $xdir] [vecscale $wy $ydir] ] 
  set corner2 [vecadd $center [vecscale $wx $xdir] [vecscale [expr (-1.0)*$wy] $ydir] ] 
  set corner3 [vecadd $center [vecscale [expr (-1.0)*$wx] $xdir] [vecscale $wy $ydir] ] 
  set corner4 [vecadd $center [vecscale [expr (-1.0)*$wx] $xdir] [vecscale [expr (-1.0)*$wy] $ydir] ] 
  set object [mol new]
  mol rename $object "$name"
  set mat [material add copy Opaque]
  material rename $mat mat$object 
  graphics $object color $color
  graphics $object material mat$object
  graphics $object triangle $corner1 $corner2 $corner3
  graphics $object triangle $corner3 $corner2 $corner4
  mol off $object
 }
 
 proc get_rot_center1 {ref_mol norv res0 res1 res2 distances} {
 
  variable pi
  # find the position of the axis, positioned in the rotation plane
  set c0 [get_xyz $ref_mol $res0]
  set c1 [get_xyz $ref_mol $res1]
  set c2 [get_xyz $ref_mol $res2]
  set v10 [vecsub $c0 $c1]
  set v12 [vecsub $c2 $c1]
  set d10 [veclength $v10]
  set d12 [veclength $v12]
  lassign $distances d0 d1 d2
 
  # get two possible positions
  # angle formed by c0, c1 and the point of interest
  set angle01p [cos_law $d10 $d1 $d0]
  #find the direction to the point, from c1
  set unit10 [vecnorm $v10]
  # rotate the unit vector +/- the angle formed.  scale by the distance from c1 and add to c1
  set p1a [vecadd $c1 [vecscale $d1 [vectrans [transabout $norv  $angle01p rad] $unit10]]]
  set p1b [vecadd $c1 [vecscale $d1 [vectrans [transabout $norv -$angle01p rad] $unit10]]]
 
  # look at next angle to find two candidates
  # angle formed by c2, c1 and the point of interest
  set angle21p [cos_law $d12 $d1 $d2]
  #find the direction to the point, from c1
  set unit12 [vecnorm $v12]
  # rotate the unit vector +/- the angle formed.  scale by the distance from c1 and add to c1
  set p2a [vecadd $c1 [vecscale $d1 [vectrans [transabout $norv  $angle21p rad] $unit12]]]
  set p2b [vecadd $c1 [vecscale $d1 [vectrans [transabout $norv -$angle21p rad] $unit12]]]
 
  # see if the two methods yielded a common point
  if {[veclength [vecsub $p1a $p2a] ] < 0.01 } {
   return $p1a
  }
 
  if {[veclength [vecsub $p1a $p2b] ] < 0.01 } {
   return $p1a
  }
 
  if {[veclength [vecsub $p2a $p1b] ] < 0.01 } {
   return $p2a
  }
 
  error "Unable to triangulate point.  Either an issue with head_d, body_d, body\[0-2\], head\[0-2\], or the reference structure.  This issue should only be possible when editing the script. If you are not editing the script, then please communicate the issue to the Whitford group."
 }
 
 proc cos_law {a b c} {
 
  # law of cosines
  set arg [expr ($a*$a+$b*$b-$c*$c)/(2*$a*$b) ]
  if {($arg > 1) || ($arg <-1)} {
   puts "cos_law issue. arg equal to $arg. This almost always arises from the values of head_d or body_d. About to fail...."
  }
  set angle [expr acos($arg)]
  return $angle
 } 
 
 proc calc_euler {ref_norv norv zero_tilt_dir ref_bond bond HEADBODY} {
  variable pi
  set proj [vecdot $ref_norv $zero_tilt_dir] 
  set zero_dir [vecnorm [vecsub $zero_tilt_dir [vecscale $proj $ref_norv]]]
  # define the basis set for expressing translations
  set transaxes {}
  lappend transaxes $zero_dir 
  lappend transaxes [vecnorm [veccross $ref_norv $zero_dir]] 
  lappend transaxes $ref_norv
  if {[veclength [veccross $ref_norv $norv ]] > 0} {
   # get the tilt direction that we will define as our zero
 
   ####### Get the line of nodes of the rotated ("norv") and reference ("ref_norv") planes
   set line_of_nodes [vecnorm [veccross $ref_norv $norv]]
   # calculate the direction of the line of nodes (tilt axis), based on our defined zero value 
   set nd [expr acos([vecdot $zero_dir $line_of_nodes])]
   set nd [expr $nd * 180 / $pi]
   set cno [veccross $zero_dir $line_of_nodes]
   if {[vecdot $cno $ref_norv] < 0} { 
    set nd -$nd
   }
   
   set phi [expr acos([vecdot $ref_bond $line_of_nodes])]
   set phi [expr $phi * 180 / $pi]
   set cno [veccross $ref_bond $line_of_nodes]
   if {[vecdot $cno $ref_norv] < 0} { 
    set phi -$phi 
   }
 
   set psi [expr acos([vecdot $line_of_nodes $bond])]
   set psi [expr $psi * 180 / $pi]
   set cno [veccross $line_of_nodes $bond]
   if {[vecdot $cno $norv] < 0} { 
     set psi -$psi 
   } 
   set theta [expr acos([vecdot $ref_norv $norv])]
   set theta [expr $theta * 180 / $pi]
 
   set rotation [expr $phi + $psi]
   set tilt $theta
   set tilt_direction $nd
   set tilt_vector $line_of_nodes
  } else {
   # This line to avoid 1 from rounding to something larger than 1 and giving an error with acos. It can really only happen for our reference configuration, but might as well take care of it.
   set pp [vecdot $ref_bond $bond]
   if { $pp > 1 } {
    set pp 1
   }
   set phi [expr acos($pp)* 180 / $pi ]; 
   set cno [veccross $ref_bond $bond]
   if {[vecdot $cno $ref_norv] < 0} { set phi -$phi }
   
   set rotation $phi
   set tilt 0
   set psi 0
   set tilt_direction 0
   set tilt_vector "not defined"
  
  }
  # only return angles from -180 to 180
  if { $rotation < -180 } { 
   set rotation [ expr $rotation+360 ]  
  }
  if { $rotation > 180 } { 
   set rotation [ expr $rotation-360 ] 
  }
 
  return [list $rotation $phi $psi $tilt $tilt_direction $tilt_vector $zero_dir $transaxes]
 }
 
 proc inter {v1 v1p v2 v2p v3} {
  # find the direction of the line of intersection of two planes.  Each set of vectors defines a plane.  Since this can produce two antiparallel solutions, take the one that has a positive projection with v3
  set dv1 [vecnorm [vecsub $v1p $v1]]
  set dv2 [vecnorm [vecsub $v2p $v2]]
  set int [vecnorm [veccross $dv1 $dv2]]
  if {[vecdot $int $v3] < 0 } {
   set int [vecscale -1 $int]
  }
  return $int
 }
 
 proc draw_arrow_scaled {mol start vec length scale color} {
  set end [vecadd $start [vecscale $length [vecnorm $vec]]]
  set merge [vecadd $start [vecscale 0.7 [vecsub $end $start]]]
  graphics $mol color $color 
  graphics $mol cylinder $start $merge radius [expr 0.04*$scale] 
  graphics $mol cone $merge $end radius [expr 0.08*$scale] 
 } 
 
 proc cleanup { string } {
  regsub -all {\}}     $string "" string
  regsub -all {\{}    $string "" string
  regsub -all {^\s+}    $string "" string
  regsub -all {\s+$}    $string "" string
  regsub -all {\s+}    $string " " string
  set a [lindex $string 0]
  set b [lindex $string 1]
  return [list $a $b]
 }
 
 proc cleanup_single { string } {
  regsub -all {\}}     $string "" string
  regsub -all {\{}    $string "" string
  regsub -all {^\s+}    $string "" string
  regsub -all {\s+$}    $string "" string
  regsub -all {\s+}    $string " " string
  return $string
 }
 
 proc format_list {list format} {
  set new ""
  foreach j $list {
   lappend new [format $format $j]
  }
  return $new
 }
 
proc print_results { optin ID_LARGE ID_SMALL CHAIN_ID_LrRNA CHAIN_ID_SrRNA body_rotation body_phi body_psi body_tilt body_tilt_direction body_tilt_vector body_translation body_origin head_rotation head_phi head_psi head_tilt head_tilt_direction head_tilt_vector head_translation head_origin RMSD_LrRNA RMSD_body RMSD_head NUM_LrRNA MATCH_LrRNA NUM_BODY MATCH_BODY NUM_HEAD MATCH_HEAD vec_body_ref vec_body_fit vec_head_ref vec_head_fit ERvec_body ERangle_body ERvec_head ERangle_head nalignLrRNA nalignSrRNA nchangeLrRNA nchangeSrRNA COR_body minimum_translation_body COR_head minimum_translation_head deltaabs_body deltaabs_head Mlarge Mbody Mhead} {
 variable radenv

 upvar $optin options
 if { [info exists radenv(ROTTEST)] } {
  set tprec 5
 } else {
  set tprec 3
 }

 set prec $options(PRECISION)
 set prec "%.${prec}f"
 if { $nalignLrRNA eq "N/A" || ![info exists options(HEADFIRST)] } {
  set head_b "N/A"
  set head_e "N/A"
 } else {
  set head_b $options(HEADFIRST)
  set head_e $options(HEADLAST)
 }
 set SSU_pdb [molinfo $ID_SMALL get filename]

 # format all lists before writing
 set head_translation [format_list $head_translation "%.${tprec}f"]
 set deltaabs_head [format_list $deltaabs_head "%.${tprec}f"]

 set head_origin [format_list $head_origin "%.8f"]
 if { $COR_head  ne "not defined" } {
  set COR_head [format_list $COR_head "%.8f"]
 }
 if {$minimum_translation_head ne "not defined" } {
  set minimum_translation_head [format_list $minimum_translation_head "%.${tprec}f"]
 }
 if {$head_tilt_vector ne "not defined"} {
  set head_tilt_vector [format_list $head_tilt_vector "%.8f"]
 }

 set head_rotation [format $prec $head_rotation]
 set head_phi [format $prec $head_phi]
 set head_psi [format $prec $head_psi]
 set head_tilt [format $prec $head_tilt]
 set head_tilt_direction [format $prec $head_tilt_direction]
 set ERangle_head [format $prec $ERangle_head]
 set vectordump ""
 set bodyvecdump ""
 set headvecdump ""
 set bodyERvecdump ""
 set headERvecdump ""
 set transMatricesdump ""

 if { [info exists options(ONLYSSU)] } {
  set lalignn ""
  set lalignrenum ""
  set matchLrRNA ""
  set matchLrRNA2 ""
  set rmsdL ""
  set ERBODY ""
  set EULBODY ""
  set BODYROT ""
  set BODYTRANS ""
  set BODY_minimum_translation_COR ""
  set BODY_minimum_translation ""
  set Lin ""
  ### extra stuff to print with -dump
  if {[info exists options(DUMPDATA)] } {
   set vectordump "    Rigid-body axes (ref model frame of ref):
        
        head (position and 3 vectors):
            reference:
            $vec_head_ref	
            model:
            $vec_head_fit	

    Translation Coordinates (internal, global values):
            head: $deltaabs_head
"
   set headvecdump "        head tilt axis      = $head_tilt_vector	
        head rot/tilt center= $head_origin\n"

   set headERvecdump "        head axis vector     = $ERvec_head 
        head axis point      = $COR_head\n"

   set transMatricesdump "
    Transformation matrices: \[ Mbody(model_body)=ref_body; Mhead(ref_head)=model_head \] 
        Mbody:  $Mbody\n
        Mhead:  $Mhead\n"

  }
 } else {

  set LSU_pdb [molinfo $ID_LARGE get filename]
  set Lin "    large rRNA : chain $CHAIN_ID_LrRNA of file $LSU_pdb (white)\n"
  set body_translation [format_list $body_translation "%.${tprec}f"]
  set deltaabs_body [format_list $deltaabs_body "%.${tprec}f"]
  set body_origin [format_list $body_origin "%.8f"]
  if { $COR_body ne "not defined" } {
   set COR_body [format_list $COR_body "%.8f"]
  }
  if {$minimum_translation_body ne "not defined" } {
   set minimum_translation_body [format_list $minimum_translation_body "%.${tprec}f"]
  }
  set body_rotation [format $prec $body_rotation]
  set body_phi [format $prec $body_phi]
  set body_psi [format $prec $body_psi]
  set body_tilt [format $prec $body_tilt]
  set body_tilt_direction [format $prec $body_tilt_direction]
  set ERangle_body [format $prec $ERangle_body]

  if {$body_tilt_vector ne "not defined"} {
   set body_tilt_vector [format_list $body_tilt_vector "%.8f"]
  }

  ### extra stuff to print with -dump
  if {[info exists options(DUMPDATA)] } {
   set vectordump "    Rigid body axes:
        body (position and 3 vectors):
            reference:
            $vec_body_ref	
            model:
            $vec_body_fit	
        
        head (position and 3 vectors):
            reference:
            $vec_head_ref	
            model:
            $vec_head_fit	

    Translation Coordinates (internal, global values):
            body: $deltaabs_body
            head: $deltaabs_head
" 
   set bodyvecdump "        body tilt axis       = $body_tilt_vector
        body rot/tilt center = $body_origin\n"

   set headvecdump "        head tilt axis       = $head_tilt_vector	
        head rot/tilt center = $head_origin\n"

   set bodyERvecdump "        body axis vector     = $ERvec_body 
        body axis point      = $COR_body\n"
   set headERvecdump "        head axis vector     = $ERvec_head 
        head axis point      = $COR_head\n"
    
   set transMatricesdump "
    Transformation matrices: \[ Mlarge(model_large)=ref_large; Mbody(ref_body)=Model_body; Mhead(ref_head)=Model_head \]
        Mlarge: $Mlarge\n
        Mbody:  $Mbody\n
        Mhead:  $Mhead\n"
  } 
  ## end extra dump stuff

  set lalignn "            large rRNA  : $nalignLrRNA\n"
  set lalignrenum "            large rRNA  : $nchangeLrRNA\n"
  set matchLrRNA "        large rRNA (blue) : $MATCH_LrRNA\n" 
  set matchLrRNA2 "        large rRNA (blue) : $MATCH_LrRNA of $NUM_LrRNA\n" 
  set rmsdL "        large rRNA = $RMSD_LrRNA\n"
  set ERBODY "        body rotation        = $ERangle_body 
$bodyERvecdump        body translation     = $minimum_translation_body\n\n"
  set EULBODY "        body: $body_phi, $body_tilt, $body_psi\n"
  set BODYROT "        body rotation        = $body_rotation
        body tilt            = $body_tilt
        body tilt direction  = $body_tilt_direction
$bodyvecdump        body translation     = $body_translation\n\n"
  set BODY_minimum_translation_COR "                body translation COR = $COR_body"
  set BODY_minimum_translation "                body translation     = $minimum_translation_body\n\n"
 }

 if { [info exists radenv(ROTTEST)] } {
  foreach name {head_tilt_vector head_origin ERvec_head COR_head body_tilt_vector body_origin ERvec_body COR_body} {
   set j [set "$name"]
   if {$j eq "not defined"} {
    set $name {0 0 0}
   }
  }
  if { [info exists options(ONLYSSU)] } {
   return [list $head_rotation $head_tilt $head_tilt_direction $head_tilt_vector $head_origin $head_translation $head_phi $head_tilt $head_psi $ERangle_head $ERvec_head $minimum_translation_head $COR_head]
  } else {
   return [list $body_rotation $body_tilt $body_tilt_direction $body_tilt_vector $body_origin $body_translation $body_phi $body_tilt $body_psi $ERangle_body $ERvec_body $minimum_translation_body $COR_body $head_rotation $head_tilt $head_tilt_direction $head_tilt_vector $head_origin $head_translation $head_phi $head_tilt $head_psi $ERangle_head $ERvec_head $minimum_translation_head $COR_head]
  } 
 } 

 if { [info exists options(PRUNE)]} {
  set pruneinfo " (after pruning, cutoff $options(PRUNEBY) Angstroms)"
 } else {
  set pruneinfo ""
 } 
 if {[info exists options(STAMP)] } {
  set aligninfo "    Alignment Information 
        Numbers of corresponding residues found via STAMP alignment
$lalignn            small rRNA  : $nalignSrRNA 
        Number of residues that were renumbered (i.e. given reference E. coli numbering)
$lalignrenum            small rRNA  : $nchangeSrRNA
"
 } else {
  set aligninfo ""
 }

 if { [info exists options(INCLUDEALL)]} { 
  set matching "    Number of residues$pruneinfo used for rigid-body alignment steps
$matchLrRNA        small body (red)  : $MATCH_BODY 
        small head (pink) : $MATCH_HEAD" 

 } else {
  set matching "    Number of matching CORE residues$pruneinfo
$matchLrRNA2        small body (red)  : $MATCH_BODY of $NUM_BODY 
        small head (pink)  : $MATCH_HEAD of $NUM_HEAD" 
 }

 puts $options(OUTPUT) "
****************** ANALYSIS SUMMARY ******************
 
INPUT MODEL	
$Lin    small rRNA : chain $CHAIN_ID_SrRNA of file $SSU_pdb (cyan)
    
CALCULATED QUANTITIES 
    Angles in degrees
    Translations given in Angstroms (internal coordinate system)
    Important Note: If you use -dump (not common), all additional 
          vectors will be in the absolute coordinate system.
$aligninfo 
$matching
        
    RMSD between model and reference configurations
$rmsdL        small body = $RMSD_body
        small head = $RMSD_head

$vectordump    
    Rotations and Translations:
      Euler Angle Based Decomposition
$BODYROT        head rotation        = $head_rotation
        head tilt            = $head_tilt
        head tilt direction  = $head_tilt_direction
$headvecdump        head translation     = $head_translation

      Euler-Rodrigues Based Decomposition (Origin - 4V9D:DA,BA): 
$ERBODY        head rotation        = $ERangle_head 
$headERvecdump        head translation     = $minimum_translation_head
$transMatricesdump
"
}
 
 proc findbesttrio {IDs COMMON delta mindist dirset} {
  # this routine searches for sets of 3 residues that can be used to define a plane, where the normal
  # vector is parallel in the rotated and unrotated conformations
  # IDs is a list of mol IDs for the unrotated/rotated conformations
  # COMMON is the list of common resids, which will be searched
  # dirset is a list of 2 residues.  the vector between them atom2-atom1 is the preferred
  # direction of the cross product. By giving atoms for the body that are pointing away from the LSU
  # and atoms for the head that are away from the body, the trio returned will also point those
  # directions
 
  # get the + direction vector
  set sel [atomselect [lindex $IDs 0]  "name P and resid [lindex $dirset 0]"]
  set pos0 [join [$sel get {x y z}]] 
  $sel delete
  set sel [atomselect [lindex $IDs 0]  "name P and resid [lindex $dirset 1]"]
  set pos1 [join [$sel get {x y z}]] 
  $sel delete
  set posvec [vecsub $pos1 $pos0]
 
  set COMLENGTH [llength $COMMON]
  set anglemin 100000
  for {set i 0} {$i < $COMLENGTH} {incr i} { 
   set ni [lindex $COMMON $i]
   for {set id 0} {$id<2} {incr id} { 
    set IDn [lindex $IDs $id]
    set sel [atomselect $IDn "resid $ni and name P"]
    set posi($IDn) [join [$sel get {x y z}]] 
    $sel delete
   } 
   set sel [atomselect [lindex $IDs 0]  "name P and resid $COMMON and (not within $mindist of (resid $ni and name P)) and (resid gt $ni)"]
   set COMMONJ [$sel get resid]
   $sel delete
   for {set j 0} {$j < [llength $COMMONJ] } {incr j} { 
    set nj [lindex $COMMONJ $j]
    for {set id 0} {$id<2} {incr id} { 
     set IDn [lindex $IDs $id]
     set sel [atomselect $IDn "resid $nj and name P"]
     set posj($IDn) [join [$sel get {x y z}]] 
     $sel delete 
     set vecij($IDn) [vecsub $posi($IDn) $posj($IDn)]
    }
    set sel [atomselect [lindex $IDs 0]  "name P and resid $COMMON and (not within $mindist of (resid $nj and name P)) and (resid gt $nj)"]
    set COMMONK [$sel get resid]
    $sel delete
    for {set k 0} {$k < [llength $COMMONK] } {incr k} { 
     set nk [lindex $COMMONK $k]
     if { ($ni+$delta <= $nj) &&  ($nj+$delta <= $nk) } {
      for {set id 0} {$id<2} {incr id} { 
       set IDn [lindex $IDs $id]
 
       set sel [atomselect $IDn "resid $nk and name P"]
       set posk [join [$sel get {x y z}]]
       $sel delete 
 
       set vecik [vecsub $posi($IDn) $posk]
       set norm($id) [vecnorm [veccross $vecik $vecij($IDn)]]
      }
      set angle [expr acos([vecdot $norm(0) $norm(1)])]
      if {$angle < $anglemin} {
       if { [vecdot $posvec $norm(0)] > 0 } {
        # already pointing in the correct direction
        puts "found new lowest tilt value: $ni $nj $nk $angle"
       } else {
        # pointing downward, so flip
        puts "found new lowest tilt value: $nk $nj $ni $angle"
       }
       set anglemin $angle
      } 
     }
    }
   }
  }
 }
 
 proc orientribosome {} {
  # this is a pre-defined view of the ribosome from the back side of the SSU
  set r "{{-0.345669 -0.768485 0.538463 0} {-0.925379 0.374278 -0.0598905 0} {-0.15551 -0.518985 -0.840518 0} {0 0 0 1}}"
  set c "{{1 0 0 40.8744} {0 1 0 -81.2588} {0 0 1 -84.263} {0 0 0 1}}"
  set s "{{0.0117307 0 0 0} {0 0.0117307 0 0} {0 0 0.0117307 0} {0 0 0 1}}"
  set g "{{1 0 0 -0.04} {0 1 0 -0.12} {0 0 1 0.62} {0 0 0 1}}"
  foreach mol [molinfo list] {
   molinfo $mol set rotate_matrix $r
   molinfo $mol set center_matrix $c
   molinfo $mol set scale_matrix $s
   molinfo $mol set global_matrix $g
  }
  display depthcue on
  display cuedensity 0.100000
  display nearclip set 0.010000
  translate by 0.000000 0.000000 -1.000000
  display projection Orthographic
  color Display Background white
  color Axes Labels black
 }
 
 proc geom_center {selection} {
  variable pi
  # set the geometrical center to 0
  set gc [veczero]
  foreach coord [$selection get {x y z}] {
   # sum up the coordinates
   set gc [vecadd $gc $coord]
  }
  return [vecscale [expr 1.0 /[$selection num]] $gc]
 }
 
 proc find_ER_between_2_systems { v1 v2 v3 v4 v5 v6 } {
  variable pi
  # this procedure calculates the Euler Rodrigues angle and vector characterizing the rotation taking the axes system (v1 v2 v3) to the axes system (v4 v5 v6): [N.B.: both systems should be left-handed orthonormal systems and the order of the vectors in both cases is x y z]
  
  set V1 [lappend v1 0]
  set V2 [lappend v2 0]
  set V3 [lappend v3 0]
  set Va {0 0 0 1}
  
  set A_tr [list "$V1" "$V2" "$V3" "$Va"];   #matrix of first system (transposed)
  
  set V4 [lappend v4 0]
  set V5 [lappend v5 0]
  set V6 [lappend v6 0]
  set Vb {0 0 0 1}
  
  set B_tr [list "$V4" "$V5" "$V6" "$Vb"];  #matrix of second system (transposed)
  
  set B [transtranspose $B_tr]
  set R [transmult $B $A_tr];      #transformation matrix (rotation matrix needed to get from the first system to the second system).
  
  set R00 [lindex $R 0 0]
  set R01 [lindex $R 0 1]
  set R02 [lindex $R 0 2]
  set R10 [lindex $R 1 0]
  set R11 [lindex $R 1 1]
  set R12 [lindex $R 1 2]
  set R20 [lindex $R 2 0]
  set R21 [lindex $R 2 1]
  set R22 [lindex $R 2 2]
  
  set Trace [expr $R00+$R11+$R22]
  if { $Trace > 2.9999992 } {
   # trace==3 means angle is less than 0.5 degrees, so everything will go haywire
   set PHI 0
   set nx 0
   set ny 0
   set nz 0
  } else {
   set PHI [expr acos(($Trace-1)/2)*180/$pi];    #Euler-Rodrigues angle
   set nx [expr ($R21-$R12)/(2*sin($PHI*$pi/180))]
   set ny [expr ($R02-$R20)/(2*sin($PHI*$pi/180))]
   set nz [expr ($R10-$R01)/(2*sin($PHI*$pi/180))]
  }
  set n_hat [list $nx $ny $nz];    #Euler-Rodrigues vector (normalized)
  if {$PHI<0} {
   set PHI [expr (-1.0)*$PHI]
   set n_hat [vecscale -1 $n_hat]
  } 
  return [list $PHI $n_hat]
 }
 
 proc findlargestblock {optin MOLID CHAIN_ID WRITE} {
  variable radenv
  upvar $optin options
  set LSEL [atomselect $MOLID "chain \"$CHAIN_ID\" and $radenv(NOTSTUFF)"]
  set last -100
  set breaks {}
  set atom 0
  set seriallist [$LSEL get serial]
  $LSEL delete
  foreach {serial} $seriallist {
   if { [expr $serial-$last] != 1 } {
    lappend breaks $atom
   }
   set last $serial
   incr atom
  }
  lappend breaks $atom
  # breaks has the first and last atoms, always.  If there are any others, then there must be something non-consecutive
  if {[llength $breaks] > 2} {
   if {$WRITE} {
    puts $options(OUTPUT) "\nNote: Chain ID $CHAIN_ID includes atoms that are no sequential in the file. This typically is due to non-rRNA atoms (often ligands) being given a redundant chain ID. For subsequent analysis, only the largest contiguous block of atoms with ID $CHAIN_ID will be considered.\n" 
   }
  }
  # if not contiguous, find largest contiguous block, rename IDs for all others and give a message
  set maxlength 0
  set maxindex 0
  for {set i 0} { $i < [expr [llength $breaks]-1] } {incr i} {
   set length [expr [lindex $breaks [expr $i+1]] - [lindex $breaks $i ]]
   if { $length > $maxlength } {
    set maxlength $length
    set maxindex $i
   }
  }
  set selt [atomselect $MOLID "serial [lindex $seriallist [lindex $breaks $maxindex]]"]
  set firstblock [$selt get serial]
  $selt delete
  set selt [atomselect $MOLID "serial [lindex $seriallist [lindex $breaks [expr $maxindex+1]]-1]"]
  set lastblock [$selt get serial]
  $selt delete
  puts $options(OUTPUT) "Largest contiguous block of atoms (indices beginning with 1): $firstblock-$lastblock
Will use this for all analysis."
  set renamelist {}
  for {set i 0} { $i < [expr [llength $breaks]-1] } {incr i} {
   if { $i != $maxindex} {
    lappend renamelist [lrange $seriallist [lindex $breaks $i] [expr [lindex $breaks [expr $i +1]]-1]  ]
   }
  }
  if { [llength $renamelist] != 0 } {
   set renamelist [join $renamelist]
   set renamesel [atomselect $MOLID "serial $renamelist"] 
   $renamesel set chain "null"
   $renamesel delete
  }
 }
 
 proc restorenumbering {ID {nodel 1}} {
  variable fullfilenames
  set origfile $fullfilenames($ID)
  set newfile  [mol new "$origfile" waitfor -1] 
  set fullfilenames($newfile) $origfile
  set orig [atomselect $ID "all"]
  set new [atomselect $newfile "all"]
 
  set M [measure fit $new $orig]
  $new move $M
  if {$nodel ==1} {
   mol delete $ID
  }
  $orig delete
  $new delete
  return $newfile
 }
 
 proc find_COR_with_minimum_absolute_translation { M domain Transformation_of_P } {
  # this procedure takes as an input a transformation matrix (4*4) that contains the rotation matrix R (3*3) and the translation vector T. The transformation matrix is applied from the origin O of the VMD reference frame(that of the reference structure pdb file). The procedure outputs a point that is the center of rotation where when R is applied from it gives the smallest absolute minimum translation. The line that passes through this point parallel to the Euler-Rodrigues vector (n_hat) is the set of all points in space that give the same minimum absolute translation if any point in it is chosen as the center of rotation from where to apply R.
 
  set R00 [lindex $M 0 0]
  set R01 [lindex $M 0 1]
  set R02 [lindex $M 0 2]
  set R10 [lindex $M 1 0]
  set R11 [lindex $M 1 1]
  set R12 [lindex $M 1 2]
  set R20 [lindex $M 2 0]
  set R21 [lindex $M 2 1]
  set R22 [lindex $M 2 2]
 
  set T0 [lindex $M 0 3]
  set T1 [lindex $M 1 3]
  set T2 [lindex $M 2 3]
 
  set T [list $T0 $T1 $T2];   #Translation vector from origin O
 
  set Trace [expr $R00+$R11+$R22]
  set PHI [expr acos(($Trace-1)/2)];    #Euler-Rodrigues angle (radians)
 
  set nx [expr ($R21-$R12)/(2*sin($PHI))]
  set ny [expr ($R02-$R20)/(2*sin($PHI))]
  set nz [expr ($R10-$R01)/(2*sin($PHI))]
  set n_hat [list $nx $ny $nz];    #Euler-Rodrigues vector (normalized)
 
 # the calculation:
  set n_cross_T [veccross $n_hat $T]
  set n_dot_T [vecdot $n_hat $T]
  set scale_1st_term [expr 1/(2*tan($PHI/2))]
 
  set first_term [vecscale $scale_1st_term $n_cross_T]
  set second_term [vecscale 0.5 [vecsub $T [vecscale $n_dot_T $n_hat]]]
  set X_perpendicular [vecadd $first_term $second_term]
 
  #COR will be chosen to be the point along the line parallel to n_hat that passes through X_perpendicular and the nearest to point P. This occurs when a=n_hat.P:
  if {$domain=="BODY"} {set P {-23.88731773 94.64822214 120.99501688}}
  if {$domain=="HEAD"} {set P {-46.9912841 101.94023883 53.25907833}}
  set P [coordtrans $Transformation_of_P $P];      ## The point P is carried with the body rotation if we are calculating COR_head and minimum_translation_vector for the head.
  
  set a [vecdot $n_hat $P]
  
  set COR [vecadd [vecscale $a $n_hat] $X_perpendicular]
  set minimum_absolute_translation_vector [vecscale $n_dot_T $n_hat]
 
  return [list $COR $minimum_absolute_translation_vector]
 }
 
 proc gettransvecs {ID_REF_LrRNA transatoms} {
  for {set i 0} { $i < 4 } {incr i} {
   set sel [atomselect $ID_REF_LrRNA "name P and resid [lindex $transatoms $i]"]
   set x($i) [join [$sel get {x y z}]]
   $sel delete
  }
  set v1 [vecnorm [vecsub $x(1) $x(0)]]
  # we will keep v2
  set v2 [vecnorm [vecsub $x(3) $x(2)]]
  set v3 [vecnorm [veccross $v1 $v2]]
  set v4 [vecnorm [veccross $v2 $v3]]
  return [list $v4 $v2 $v3]
 }
 
 proc transvecstrans {translation tvecs} {
  # translation is a vector describing a translational motion
  # tvecs is the new basis
  set v0 [lindex $tvecs 0]
  set v1 [lindex $tvecs 1]
  set v2 [lindex $tvecs 2]
 
  set proj0 [vecdot $translation $v0]
  set proj1 [vecdot $translation $v1]
  set proj2 [vecdot $translation $v2]
 
  return [list $proj0 $proj1 $proj2]
 }
 
 proc drawtransvecstrans {tvecs BH} {
  variable labels
  # translation is a vector describing a translational motion
  # tvecs is the new basis
  set v0 [lindex $tvecs 0]
  set v1 [lindex $tvecs 1]
  set v2 [lindex $tvecs 2]
 
  draw_fitted_axis {0 0 0} $v0 "red" $labels(tv1$BH)  
  draw_fitted_axis {0 0 0} $v1 "white" $labels(tv2$BH) 
  draw_fitted_axis {0 0 0} $v2 "blue" $labels(tv3$BH)
 }
 
 proc closelogs {optin} {
  upvar $optin options
  if { $options(OUTPUT) ne "stdout" } {
   close $options(OUTPUT)
  }
  if { $options(ERROR) ne "stderr" } {
   close $options(ERROR)
  }
 }
 
 proc finalwrapup {optin} {
  variable radenv
  upvar $optin options
  variable RTD
  foreach i { ID_LARGE ID_SMALL CHAIN_ID_LrRNA CHAIN_ID_SrRNA ID_REF_SrRNA_5 } { 
   set $i $RTD($i)
  }
  # During alignment we internally renumber a range of things.
  # here, we are going to re-load the original structures, align them
  # and recolor them, so that querying atoms will not give unusual values
  if {$ID_LARGE == $ID_SMALL} {
   set same 0
  } else {
   set same 1
  }
  puts "Visualization Information:
    Will reload and reorient the input model"
  if {[info exists options(SAVENAME)] && ! [info exists options(ONLYLSU)] } {
   # reset the numbering and reload the PDB
   set ID_SMALLt [restorenumbering $ID_SMALL 0]
   # findlargestblock finds the largest continuous set of residues with the given chain ID
   # all others that have the same ID are set to ID null
   findlargestblock options $ID_SMALLt $CHAIN_ID_SrRNA 0
   set tmpsel [atomselect $ID_SMALLt "chain \"$CHAIN_ID_SrRNA\" "]
   $tmpsel  writepdb "$options(SAVENAME)_SSU.pdb"
   $tmpsel delete
   mol delete $ID_SMALLt
  }
  if { ! [info exists options(ONLYLSU)] } {
   set ID_SMALL [restorenumbering $ID_SMALL]
   set CURHISTORY [split [lindex $radenv(RADtoolhistory) end] "\ "] 
   lset CURHISTORY end $ID_SMALL 
   drawreptubelines $ID_SMALL $CHAIN_ID_SrRNA 1.3 21
   puts "        Molecule ID of the original SSU structure (after alignment) is $ID_SMALL"
  }
  if { ! [info exists options(ONLYSSU)] } {
   if {$same == 0} {
    puts "        Molecule ID of the LSU structure is $ID_SMALL"
    set ID_LARGE $ID_SMALL
   } else {
    set ID_LARGE [restorenumbering $ID_LARGE]
    puts "        Molecule ID of the original LSU structure (after alignment) is $ID_LARGE"
   }
   if {[info exists options(SAVENAME)]} {
    # reset the numbering and reload the PDB
    set ID_LARGEt [restorenumbering $ID_LARGE 0]
    # findlargestblock finds the largest continuous set of residues with the given chain ID
    # all others that have the same ID are set to ID null
    findlargestblock options $ID_LARGEt $CHAIN_ID_LrRNA 0
    set tmpsel [atomselect $ID_LARGEt "chain \"$CHAIN_ID_LrRNA\" "]
    $tmpsel  writepdb "$options(SAVENAME)_LSU.pdb"
    $tmpsel delete
    mol delete $ID_LARGEt
   }
  
   if { ! [info exists options(ONLYLSU)] } {
    lset CURHISTORY end-1 $ID_LARGE
   } 
   drawreptubelines $ID_LARGE $CHAIN_ID_LrRNA 1.3 8 $same
  }
  if { ! [info exists options(ONLYLSU)] } {
   lset radenv(RADtoolhistory) end $CURHISTORY 
   puts "    If you choose to use the -animate2 option, this ribosome is index [expr [llength $radenv(RADtoolhistory)]-1]"
   mol top $ID_SMALL
   if {[info exists options(ANIMATE)]} {
    # turn on animated frames
    mol on $ID_REF_SrRNA_5
    mol top $ID_REF_SrRNA_5
    mol rename $ID_REF_SrRNA_5 "animation"
   } else {
    mol delete $ID_REF_SrRNA_5
   }
  }
  orientribosome
 }

 proc description {} {
  variable radversion
  variable labels
 # this just returns a description of RADtool stuff.  It can be called anywhere we want to give info.
  puts "
******************************   Description of RADtool   ******************************


RADtool is a plugin for VMD that is specifically designed to analyze 
structures of the ribosome. The tool allows one to provide a single 
structure file, or several files, that contain the rRNA of a ribosomal
large subunit (LSU), small subunit (SSU), or both. The tool can perform
structure alignment of the LSU, SSU (head and body) and then calculate
the Euler angles associated with the orientation of the body and/or head.

As output, the script will print (or save to file), the rotation angle,
the tilt angle and the tilt direction. In addition, if the analyzed 
structure can not be related to a reference E. coli structure through pure
rigid-body rotations, then the required translation is also calculated.


                               VISUALIZATION INFORMATION
This tool will generate many forms of visualizations, which include
representations of the LSU and SSU ``cores'', rotation vectors, and more.
Not all items are shown by default. If they are not shown, double click
the red \"D\" next to the name in the VMD Main window. While all items may
be altered by the user, the default display of each is described below.
Labels are the \"Molecule\" names that are shown in the VMD Main window.

  STRUCTURAL MODELS:
    $labels(arl) - reference E. coli LSU rRNA, aligned to the input model.
        LSU core shown in blue tubes 
    $labels(aib) - reference E. coli SSU body, aligned to the input model
        SSU body core shown in red tubes
    $labels(aih) - reference E. coli SSU head, aligned to the input model
        SSU head core shown in pink tubes
    $labels(rbp) - the reference E. coli position of the SSU body, relative to
        the position of the LSU in the input model SSU body shown in yellow tubes
    $labels(rhp) - the reference E. coli position of the SSU head, relative to
        the position of the SSU body in the input model SSU head shown in orange tubes
    $labels(ani) - animation showing body rotation, tilt and translation,
        as well as head rotation, tilt and translation. Animation begins from the 
        reference E. coli model. This is only generate with the -animate option 
        (command-line option). When using the GUI, select \"Animate from reference\" 

  ROTATION AXES AND PLANES: 
   Note: Displayed axes are centered at the point about which translation
        is minimized 
    $labels(brar) - primary body rotation axis, fixed to reference model  
    $labels(bram) - primary body rotation axis, fixed to input model
    $labels(brpr) - primary body rotation plane, fixed to reference model
        This plane is defined to be perpendicular to $labels(brar) 
    $labels(brpm) - primary body rotation plane, fixed to input model
        This plane is defined to be perpendicular to $labels(bram) 
    $labels(bzta) - defined tilt direction of zero for the body
    $labels(bta)  - calculated body tilt axis. 
        This is parallel to the cross product of the primary body rotation
        axes of the reference and input models. It is also defined as the
        intersection of $labels(brpr) and $labels(brpm)

    $labels(hrar) - primary head rotation axis, fixed to reference model  
    $labels(hram) - primary head rotation axis, fixed to input model
    $labels(hrpr) - primary head rotation plane, fixed to reference model
        This plane is defined to be perpendicular to $labels(hrar) 
    $labels(hrpm) - primary head rotation plane, fixed to input model
        This plane is defined to be perpendicular to $labels(hram) 
    $labels(hzta) - defined tilt direction of zero for the head
    $labels(hta)  - calculated head tilt axis. 
        This is parallel to the cross product of the primary head rotation
        axes of the reference and input models. It is also defined as the
        intersection of $labels(hrpr) and $labels(hrpm)
 
    $labels(erbra) - Euler-Rodrigues rotation axis, body
    $labels(erhra) - Euler-Rodrigues rotation axis, head

  TRANSLATION AXES:
   Translations are given in an internally-defined coordinate system. Note
        that the translation vectors are drawn to scale, which means they
        are typically only a few Angstroms long.
    $labels(btran) - body translation vector
    $labels(htran) - head translation vector

   The directions of the three basis vectors  (v1, v2, 3) can be shown.
   
   body translation basis vectors  
    $labels(tv1body) - vector 1 (\"X\" axis)
    $labels(tv2body) - vector 2 (\"Y\" axis)
    $labels(tv3body) - vector 3 (\"Z\" axis)

   head translation basis vectors  
    $labels(tv1head) - vector 1 (\"X\" axis)
    $labels(tv2head) - vector 2 (\"Y\" axis)
    $labels(tv3head) - vector 3 (\"Z\" axis)


                                  ANALYSIS INFORMATION

RADtool will calculate various quantities, which are either printed to the 
Tk Console, or written to file. Here are some descriptions of each quantity.

  Rotations and Translations:
   Euler Angle Based Decomposition
    body rotation        - rotation about $labels(brar)
    body tilt            - tilt rotation about $labels(bta)
    body tilt direction  - the direction of $labels(bta)
    body translation     - any orientation change that can not be
        accounted for by rotation

    head rotation        - rotation about $labels(hrar)
    head tilt            - tilt rotation about $labels(hta)
    head tilt direction  - the direction of $labels(hta)
    head translation     - any orientation change that can not be
        accounted for by rotation

   Euler-Rodrigues Based Decomposition (Origin - 4V9D:DA,BA): 
    body rotation        - rotation about $labels(erbra) 
    body translation     - required translation

    head rotation        - rotation about $labels(erhra)
    head translation     - required translation


                                SOME TIPS

- Any molecule shown can be saved to file by right clicking on the name
  in the VMD Main window, and selecting \"Save Coordinates\"

- The default representation of the input model is to only show the rRNA.
  However, the full input structure is in memory.  To display the full
  input model, go to the VMD Main->Graphics->Representations window,
  select the input model from the drop-down menu and set the \"Selected Atoms\" 
  to \"all\"


For any questions about RADtool, please contact:
Whitford Research Group
Center for Theoretical Biological Physics
Northeastern University
p.whitford@northeastern.edu
"
 }
}

set ::RADTOOL::radversion "1.2beta"
set ::RADTOOL::radenv(ROTATIONPATH) [ file dirname [ file normalize [ info script ] ] ]
::RADTOOL::setrefslabels
source $::RADTOOL::radenv(ROTATIONPATH)/single_align.tcl
source $::RADTOOL::radenv(ROTATIONPATH)/RAD-animate.tcl
source $::RADTOOL::radenv(ROTATIONPATH)/RAD-bundle.tcl
source $::RADTOOL::radenv(ROTATIONPATH)/RAD-tests.tcl
 
 
