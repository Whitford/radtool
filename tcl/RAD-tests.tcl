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


namespace eval RTTESTS {
 ################# Testing routines #################
 
 proc RT-tests args {
  set thresh 0.001
  set trajthresh 0.02
  if { $args eq "all" } {
   puts stderr "Will runs all tests, including trajectories (may be slow)."
  } elseif {  $args eq "trajonly" } {
   puts stderr "Will runs all trajectory tests, only (may be slow)."
  } elseif {  $args eq "" } {
   puts stderr "Will runs tests for single structures (not testing trajectories)."
  } else {
   puts stderr "\"$args\" is not a supported option with RADtool -test. May run without an option, or with \"all\" or \"trajonly\". These options are only necessary if you plan to analyze trajectories."
   return 
  } 
  set SCRIPTPATH $::RADTOOL::radenv(ROTATIONPATH)
  set origDir [pwd]
  set tdir "$SCRIPTPATH/../testingdir"
  set PDBDIR "$SCRIPTPATH/../share/testfiles"
  set ::RADTOOL::radenv(ROTTEST) 1
  set passed {} 
  set failed {}
  set numtests 0

  if { ! [file exists $tdir]} {
   file mkdir $tdir
  } 
  cd $tdir
 
  if { $args ne "trajonly" } {
   set tests [open "$PDBDIR/testlist" r ]
   set data [read $tests]
   close $tests
   set testlist [split $data \n]
   foreach {test} $testlist {
    regsub -all {^\s+}  $test {} test
    regsub -all {\s+$}  $test {} test
    regsub -all {\s+\|}  $test {|} test
    regsub -all {\|\s+}  $test {|} test
    regsub -all {\s+}  $test { } test
 
    set thistest [split $test "\|"]
 
    if { [llength $thistest] == 0} {
     continue
    }
 
    set name [lindex $thistest 0]
 
    if { [llength $thistest] > 6 } {
     puts stderr "Error: testlist file has line with too many entries: $thistest"
     cd $origDir
     return
    }
 
 
    set matchhead [regexp {(-head$)} $name]
 
    if { $matchhead == 0 } {
     incr numtests
     if { [llength $thistest] < 6 && $name ne "bundle" && $name ne "bundlecif" } {
      puts "warning: testlist file has line with too few entries: $thistest"
      puts "skipping test"
      continue
     }
     if { $name eq "bundle" || $name eq "bundlecif"  } {
      set structures [lindex $thistest 1]
      set filelist {}
      foreach {test} $structures {
       lappend filelist "$PDBDIR/$test " 
      }
      set filelist [join $filelist]
      set vals [lindex $thistest 2]
     } else {
      set l [lindex $thistest 1]
      set s [lindex $thistest 2]
      set lc [lindex $thistest 3]
      set sc [lindex $thistest 4]
      set vals [lindex $thistest 5]
     } 
     puts "Starting test: $name"
     if { $name eq "basic" } {
      set values [callRAD -precision 2 -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -noprune -notall -overwrite]
     } elseif { $name eq "bundle" || $name eq "bundlecif" } {
      set values [callRAD  -precision 2 -stamp -bundle $filelist -o $tdir/ribrottesting.tmpfile -overwrite]
     } elseif { $name eq "alignin" } {
      set values [callRAD  -precision 2 -align_in $PDBDIR/$l.align -noprune -notall  -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } elseif { $name eq "stamp" || $name eq "stamp-single" } {
      set values [callRAD  -precision 2 -nofindhead -stamp -noprune -notall -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite -align_out $tdir/ao.out -cores_out $tdir/co.out]
      catch {file delete -force $tdir/ao.out} tv
      catch {file delete -force $tdir/co.out } tv
     } elseif { $name eq "findhead"} {
      set values [callRAD  -precision 2 -stamp -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite -align_out $tdir/ao.out -cores_out $tdir/co.out]
      catch {file delete -force $tdir/ao.out} tv 
      catch {file delete -force $tdir/co.out} tv 
     } elseif { $name eq "findhead-outin"} {
      set values [callRAD  -precision 2 -stamp -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -cores_out $tdir/tmp.co -align_out $tdir/tmp.ao -overwrite]
      set values [callRAD  -precision 2 -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -cores_in $tdir/tmp.co -align_in $tdir/tmp.ao -overwrite]
      removetmpfiles $tdir
     } elseif { $name eq "findhead-outin-sub"} {
      # this test will perform alignment for the LSU-SSU assembly, but then repeat only the SSU analysis using the alignment and cores files from the full assembly.
      set values [callRAD  -precision 2 -stamp -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -cores_out $tdir/tmp.co -align_out $tdir/tmp.ao -overwrite]
      set values [callRAD  -precision 2 -s $PDBDIR/$s -SsUonly -sc $sc -o $tdir/ribrottesting.tmpfile -cores_in $tdir/tmp.co -align_in $tdir/tmp.ao -overwrite]
      removetmpfiles $tdir
 
     } elseif { $name eq "alignin-prune" } {
      set values [callRAD  -precision 2 -align_in $PDBDIR/$l.align -pruneby 1 -notall -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } elseif { $name eq "alignin-prune-all"} {
      set values [callRAD  -precision 2 -align_in $PDBDIR/$l.align -pruneby 1 -l $PDBDIR/$l -s $PDBDIR/$s -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } else {
      puts "ERROR: name $name not recognized"
      cd $origDir
      return
     }
 
     lassign [compareline $vals $values $thresh] pass difflist
 
     if {$pass == 1} {
      puts "Passed test: $name\n"
      lappend passed $name
     } else {
      puts stderr "Failed test: $name\nexpected:\n$vals\nobtained:\n$values"
      if { $values ne {} } {
       puts stderr "The following elements differed from the benchmark by more than the allowable threshold ($thresh for quanties smaller than 50, or [expr $thresh*10] otherwise):\n $difflist\n\n"
      }
      lappend failed $name
     }
     removertmpfile $tdir
 
    } elseif { $matchhead == 1 } {
     incr numtests
     if { [llength $thistest] != 4 } {
      puts "warning: testlist file has wrong number of entries: $thistest"
      puts "skipping test"
      continue
     }
     set s [lindex $thistest 1]
     set sc [lindex $thistest 2]
     set vals [lindex $thistest 3]
     puts "Starting test: $name"
 
     if { $name eq "basic-head" } {
      set values [callRAD  -precision 2 -SSUonly -s $PDBDIR/$s  -sc $sc -o $tdir/ribrottesting.tmpfile -noprune -notall -overwrite]
     } elseif { $name eq "alignin-head" } {
      set values [callRAD  -precision 2 -SSUonly -align_in $PDBDIR/$s.align -noprune -notall  -s $PDBDIR/$s -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } elseif { $name eq "stamp-head" } {
      set values [callRAD  -precision 2 -nofindhead -stamp  -s $PDBDIR/$s -SSUonly -noprune -notall -sc $sc -o $tdir/ribrottesting.tmpfile -cores_out $tdir/co.out -align_out $tdir/ao.out -overwrite]
     catch {file delete -force $tdir/ao.out} tv 
     catch {file delete -force $tdir/co.out} tv 
     } elseif { $name eq "findhead-outin-head"} {
      set values [callRAD  -precision 2 -stamp -s $PDBDIR/$s -ssuonly -sc $sc -o $tdir/ribrottesting.tmpfile -cores_out $tdir/tmp.co -align_out $tdir/tmp.ao -overwrite]
      set values [callRAD  -precision 2 -s $PDBDIR/$s -SsUonly -sc $sc -o $tdir/ribrottesting.tmpfile -cores_in $tdir/tmp.co -align_in $tdir/tmp.ao -overwrite]
      removetmpfiles $tdir
     } elseif { $name eq "stamp-outin-head" } {
      set values [callRAD  -precision 2 -nofindhead -stamp  -s $PDBDIR/$s -SSUonly -noprune -notall -sc $sc -o $tdir/ribrottesting.tmpfile -align_out $tdir/tmp.ao -cores_out $tdir/tmp.co -overwrite]
      set values [callRAD  -precision 2 -s $PDBDIR/$s -SSUonly -noprune -notall -sc $sc -o $tdir/ribrottesting.tmpfile -align_in $tdir/tmp.ao -cores_in $tdir/tmp.co -overwrite]
      removetmpfiles $tdir
     } elseif { $name eq "alignin-prune-head" } {
      set values [callRAD  -precision 2 -align_in $PDBDIR/$s.align -pruneby 1 -SSUonly -notall -s $PDBDIR/$s -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } elseif { $name eq "alignin-prune-all-head"} {
      set values [callRAD  -precision 2 -align_in $PDBDIR/$s.align -pruneby 1 -SSUonly -s $PDBDIR/$s -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite]
     } else {
      puts "ERROR: name $name not recognized"
      cd $origDir
      return
     }
 
     lassign [compareline $vals $values $thresh] pass difflist
 
     if {$pass == 1} {
      puts "Passed test: $name\n"
      lappend passed $name
     } else {
      puts stderr "Failed test: $name\nexpected:\n$vals\nobtained:\n$values"
      puts stderr "The following elements differed from the benchmark by more than the allowable threshold ($thresh for quanties smaller than 50, or [expr $thresh*10] otherwise):\n [join $difflist]\n\n"
      lappend failed $name
     }
     removertmpfile $tdir
    } 
    mol delete all
   }
  }
  if { $args eq "all" || $args eq "trajonly" } {
   set ::RADTOOL::radenv(ROTTESTTRAJ) 1
   # format of the file is to list the method on a line, followed by the PDB and traj files
   # next line, the flags and number of frames (N, values given, below)
   # next N lines give the values for each frame
   # stamp|6gz4-pdb-bundle1.pdb|6gz4-pdb-bundle2.pdb|G|A|-1.56 3.60 -19.61 19.42 5.20 -55.24 
   set tests [open "$PDBDIR/testlist.traj" r ]
   set data [read $tests]
   close $tests
   set testinfo [split $data \n]
   set nlines [ llength $testinfo ] 
   set trajnum 0
   for {set i 0} {$i < $nlines-1} {incr i} {
    incr trajnum
    incr numtests
    set type [lindex $testinfo $i]
    puts "testing trajectory $type"
    set name "traj.$type"
    incr i
    if { $type eq "skip" } {
     set skip [lindex $testinfo $i]
     incr i
    } elseif { $type eq "segment" } {
     set segment [lindex $testinfo $i]
     incr i
    }
    set files [lindex $testinfo $i]
    set pdb [lindex [split $files "\ "] 0]
    set filelist ""
    for {set j 0} {$j < [llength $files]} {incr j} {
     set filelist "$filelist $PDBDIR/[lindex $files $j]"
    }
    incr i
    set chains [ lindex $testinfo $i ]
    set sc [lindex [split $chains "\ "] 0]
    set lc [lindex [split $chains "\ "] 1]
    if { $type eq "full" } {
     set values [callRAD  -precision 2 -align_in $PDBDIR/$pdb.align  -traj $filelist -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite -ot $tdir/ribrottesting.traj.out.tmpfile]
    } elseif { $type eq "skip" }  {
     set values [callRAD  -precision 2  -align_in $PDBDIR/$pdb.align -t_step $skip -traj $filelist -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite -ot $tdir/ribrottesting.traj.out.tmpfile]
    } elseif { $type eq "segment" }  {
     set values [callRAD  -precision 2  -align_in $PDBDIR/$pdb.align -t_seglength $segment -traj $filelist -lc $lc -sc $sc -o $tdir/ribrottesting.tmpfile -overwrite -ot $tdir/ribrottesting.traj.out.tmpfile]
    } else {
     puts stderr "Error:  $traj option \"$type\" not supported for testing\nExiting without completing.\n"
     cd $origDir
     return
    }
 
    incr i
    set refvals [lindex $testinfo $i]
 
    set fileh [open "$tdir/ribrottesting.traj.out.tmpfile" r ]
    set trajvals [read $fileh]
    close $fileh
 
    set fileh [open "$PDBDIR/$refvals" r ]
    set targetvals [read $fileh]
    close $fileh
 
    set passthis 0
 
    if { $targetvals eq $trajvals } {
     #output is identical
     set passthis 1
    } else {
     # output is not identical
     set datatraj [split $trajvals "\n"]
     set datatarget [split $targetvals "\n"]
     if { [llength $datatraj] == [llength $datatarget] } {
      # will check line-by-line, if they are the same length
      set sameline 0
      for {set j 1 } { $j < [llength $datatraj] } {incr j} {
       # skipping first line, since it is supposed to be a comment
       set trajline [join [lindex $datatraj $j]]
       set targetline [join [lindex $datatarget $j]]
       lassign [compareline $targetline $trajline $trajthresh] pass difflist
       if { $pass == 1 } {
        incr sameline
       } else {
        puts stderr "Difference in frame # $j, entry/entries # [join $difflist]"
       }
      }
      if { [expr $sameline + 1] == [llength $datatraj] } {
       set passthis 1
      }
     } 
    }
 
    if {$passthis == 1} {
     puts "Passed test: $name  "
     lappend passed $name
    } else {
     puts stderr "Failed test: $name
 expected:\n$targetvals
 calculated:\n$trajvals
 Note: Numerical/rounding differences that are platform-dependent may lead to
 occasional test failures.  If you are interested in using this tool for the 
 analysis of simulations, and some tests are failed, you should decide if the 
 differences between expected and calculated values are significant. For example, 
 a single quantity being off by 0.01 is likely nothing to be concerned about. 
 But, if an angle were to differ by 0.5, there is likely a problem with the tool."
     lappend failed $name
    }
    removertmpfile $tdir
    removerttmpfile $tdir
    mol delete all
   }
   unset ::RADTOOL::radenv(ROTTESTTRAJ)
  }
 
  if { [llength $passed] == $numtests } {
   if { $args eq "all" } {
    puts stderr "\n\nRADTOOL-TEST PASSED ALL TESTS!\n\n"
   } elseif { $args eq "trajonly" } {
    puts stderr "\n\nRADTOOL-TEST PASSED ALL TRAJECTORY TESTS!\n\n"
   } else {
    puts stderr "\n\nRADTOOL-TEST PASSED ALL (non-trajectory-based) TESTS!\n\n"
   }
  } else {
 
   puts stderr "RADTOOL-TEST Summary: Passed [llength $passed] tests and failed [llength $failed] tests."
  
   if { [llength $passed] > 0 } {
    puts stderr "	Passed tests: $passed"
   }
   if { [llength $failed] > 0 } {
    puts stderr "	Failed tests: $failed"
   }
  }
  unset ::RADTOOL::radenv(ROTTEST)
  cd $origDir
 }
 
 ##############END OF TESTING ROUTINES###############
 
 proc compareline {vals values thresh} {
  set pass 0
  set expected [join $vals]
  set obtained [join $values]
  set difflist {}
  if { [llength $expected] == [llength $obtained] } {
   set same 0
   for {set i 0} {$i < [llength $expected]} {incr i} {
    set v1 [lindex $expected $i]
    set v2 [lindex $obtained $i]
    if {$v1 == $v2} {
     incr same
    } elseif {[expr sqrt((($v1-$v2))**2) ] < $thresh } {
     # if not identical, make sure they are close enough
     incr same
    } elseif {[expr sqrt((($v1-$v2))**2) ] < [expr 10*$thresh] && [expr sqrt($v1**2)] > 50.0 } {
     # since it is a larger quantity, we can relax the precision of the comparison
     incr same
    } else {
     lappend difflist "\n\t[expr $i+1] ($v1 vs. $v2)"
    }
   }
   if {$same == [llength $expected]} {
    set pass 1
   }
  }
  return [list $pass $difflist]
 }
 
 proc callRAD {args} {
  set values {}
  if {[catch {
   set values [::RADTOOL::RADtool $args]
  } errorMessage] != 0} {
   puts stderr "RADtool exited with the following error:\n$errorMessage"
  }
  return $values
 }
 
 proc removetmpfiles {dir} {
   if {[catch {
    file delete -force $dir/tmp.ao 
    file delete -force $dir/tmp.co 
   } errorMessage] != 0} {
    puts stderr "unable to remove tmp files.\n$errorMessage"
   }
 }
 
 proc removertmpfile {dir} {
   if {[catch {
    file delete -force $dir/ribrottesting.tmpfile
   } errorMessage] != 0} {
    puts stderr "unable to remove tmp file.\n$errorMessage"
   }
 } 
 
 proc removerttmpfile {dir} {
   if {[catch {
    file delete -force $dir/ribrottesting.traj.out.tmpfile 
   } errorMessage] != 0} {
    puts stderr "unable to remove tmp file.\n$errorMessage"
   }
 } 

}  
