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

namespace eval RTBUNDLE {

 proc RT-bundle {args} { 
  
  splashmessageBUNDLE
  set exec "RADtool"
 
  lassign [parse_argsBUNDLE [join $args]] pdblist RADtoolflags outfiles asubunit
  if { $pdblist eq "\.help\." } {
   return
  }
  lassign [loadPDBs $pdblist] names IDs
 
  if {$asubunit == 0 } {
   set candidates [getCandidatesComplete $IDs]
   set val [runComplete ::RADTOOL::$exec $outfiles $names $candidates $RADtoolflags]
   if { [info exists ::RADTOOL::radenv(ROTTEST)] } {
    return $val
   }
  } else {
   set candidates [getCandidatesSingle $IDs $asubunit]
   runSingleSU ::RADTOOL::$exec $outfiles $names $candidates $RADtoolflags $asubunit
  }
  # clean up
  foreach {tid} $IDs {
   mol delete $tid
  }
 }
 
 ################# END OF MAIN ######################
 
 
 #Simpler and seems to work with https and all OS.
 proc mcurl { url file } {
  if {[catch {exec curl --silent --show-error --fail -o $file $url} error]} {
   puts "$file could not be retrieved from $url due to the following error:
 $error
 $file was not downloaded. 
 " 
  } else {
   puts "$file downloaded."
  }
 }
 
 proc download_bundle { pdbdownload } {
  puts "
 Will try to download $pdbdownload"
  set dir1 [string range $pdbdownload 1 2]
  set dir2 [string range $pdbdownload 0 3]
  set url [format "https://files.rcsb.org/pub/pdb/compatible/pdb_bundle/%s/%s/%s" $dir1 $dir2 $pdbdownload]
  mcurl $url $pdbdownload
  if {  [file exists $pdbdownload] > 0 } {
   exec tar xf $pdbdownload 
   set fn 1
   if { [file exists $dir2-pdb-bundle$fn.pdb] > 0 } {
    while { [file exists $dir2-pdb-bundle$fn.pdb] > 0} {
     lappend pdblist $dir2-pdb-bundle$fn.pdb
     incr fn
    }
    return $pdblist
   } else {
    puts "$pdbdownload could be downloaded for RCSB database, though it is not a valid file."
   }
  } else {
   puts "Unable to download $pdbdownload"
  }
 }
 
 proc download_single_pdb { pdbdownload_orig } {
  puts "
 Will try to download $pdbdownload_orig"
  set url [format "https://files.rcsb.org/download/%s" $pdbdownload_orig]
  mcurl $url $pdbdownload_orig
  if { [file exists $pdbdownload_orig] eq 0 } {
   puts "Unable to download $pdbdownload_orig"
  } else {
   lappend pdblist $pdbdownload_orig
   return $pdblist
  }
 }
 
 proc download_single_cif { pdbdownload_cif } {
  puts "
 Will try to download $pdbdownload_cif" 
  set url [format "https://files.rcsb.org/download/%s" $pdbdownload_cif]
  mcurl $url $pdbdownload_cif
  if {  [file exists $pdbdownload_cif] eq 0} {
   puts "Unable to download $pdbdownload_cif"
  } else {
   lappend pdblist $pdbdownload_cif
   return $pdblist
  }
 }
 
 proc getpdb {pdbdownload} {
  if { [string length $pdbdownload] != 4 } {
   error "When using the -download option, you must provide the 4-letter code for the RCSB entry."
  }
  puts "-download flag given."
  set curdir [pwd]
  if {![file writable $curdir]} {
   error "Unable to download structures to the current directory. You do not appear to have write access. Current directory: $curdir\nPerhaps try working/saving to a different directory."
   return
  }
  puts "Checking whether $pdbdownload structure has already been downloaded."
  set pdbdownload [string tolower $pdbdownload]
  set pdbdownload_orig [format "%s.pdb" $pdbdownload]
  set pdbdownload_cif [format "%s.cif" $pdbdownload]
  set pdbdownload [format "%s-pdb-bundle.tar.gz" $pdbdownload]
  if {[info exists radenv(GUION)]} {
   set STR " using Option 2"
  } else {
   set STR " using the -bundle flag"
  }
  if { [file exists $pdbdownload] > 0 } {
   error "$pdbdownload is already present in the current directory. This is probably the desired structure, since PDB bundle names are quite specific. However, just to be safe, the script is going to exit. You can always analyze structures that you have already downloaded$STR."
  } elseif { [file exists [string trim $pdbdownload ".gz"]] > 0 } {
   error "[string trim $pdbdownload ".gz"] is already present in the current directory. This is probably the desired structure, since PDB bundle names are quite specific. However, just to be safe, the script is going to exit. You can always analyze structures that you have already downloaded$STR."
  } elseif { [file exists $pdbdownload_cif] > 0 } {
   error "$pdbdownload_cif is already present in the current directory. This is probably the desired structure, since CIF names are quite specific. However, just to be safe, the script is going to exit. You can always analyze structures that you have already downloaded using$STR."
  } else {
   set choice "1"
   while { ![info exists pdblist] || $pdblist eq "" && $choice < 6 } {
    if { $choice eq 1 || $choice eq 4 } {
     set pdblist [download_bundle $pdbdownload]
    }
    if { $choice eq 2 || $choice eq 5 } {
     set pdblist [download_single_pdb $pdbdownload_orig]
    }
    if { $choice eq 3 } {
     set pdblist [download_single_cif $pdbdownload_cif]
    }
    incr choice
   }
   if { [file exists $pdbdownload] > 0 && $pdblist eq "" } {
    error "$pdbdownload was downloaded, but it doesn't appear to contain readable structural information. This sometimes happens when a file is corrupted on the RCSB database. We suggest trying to manually download and see if the structure is readable."
  } elseif { [file exists [string trim $pdbdownload ".gz"]] > 0  && $pdblist eq "" } {
    error "[string trim $pdbdownload ".gz"] is in the current directory, but it doesn't appear to contain readable structural information. This sometimes happens when a file is corrupted on the RCSB database. We suggest trying to manually download and see if the structure is readable."
  } elseif { [file exists $pdbdownload_cif] > 0  && $pdblist eq "" } {
    error "$pdbdownload_cif is in the current directory, but it doesn't appear to contain readable structural information. This sometimes happens when a file is corrupted on the RCSB database. We suggest trying to manually download and see if the structure is readable."
  } elseif { $pdblist eq "" } {
    error "
 
 Downloading the structure file from RSCB has failed. This could happen for a variety of reasons, including:
     - You do not have write access to the working directory
     - The file in the RCSB database is corrupted
     - The accession code ([string range $pdbdownload 0 3]) is not valid 
     - A connection issue, such as a firewall, has blocked VMD from connecting to RCSB
     - RCSB is temporarily not accessible
 
 It may be necessary to manually download the structure file and then provide it to RADtool via the command line, or GUI.
 
 "
    }
  }
  return $pdblist
 }
 
 proc splashmessageBUNDLE {} {
 
  puts "
 
 Ribosome Angle Decomposition Bundle Tool being called  
 
 Will take structure files for a ribosome, and then try to 
 identify the LSU and/or SSU rRNA. Once identified, will try to
 calculate the rotation, tilt, tilt directon and translations,
 unless -LSUonly is given.
 
 "
 }
 
 proc parse_argsBUNDLE {args} {
 # flatten the arg list
  set args [join $args]
  set argl -1
  while { $argl != [llength $args]} {
   set argl [llength $args]
   set args [join $args]
  }
  set asubunit 0
  set supopts [ list "bundle" "download" "e" "o" "oe" "ot" "align_out" "cores_out" "overwrite" "dump" "stamp" "animate" "notall" "noprune" "findhead" "pruneby" "ssuonly" "lsuonly" "precision" "findhead2"]
  foreach {opt} $supopts {
   set options(-$opt) 0
  }
  set singlargs [ list "e" "o" "oe" "ot" "align_out" "cores_out" "pruneby" "precision" "animate"]
  foreach {opt} $singlargs {
   set singa(-$opt) 0
  }
  
  set saveflag 0
  set saveval 0
  set downloadflag 0
  set RADtoolflags {}
  set pdblist {}
  set pdbdownload ""
  foreach {nameo} $args {
   if { [regexp -nocase {(^-bundle$)} $nameo] == 1 } {
    set bundleon 0
   }
   set name [string tolower $nameo]
   if { [regexp {(^-)} $nameo] == 1 && ! [info exists options($name)] } {
    error "Flag not supported with -bundle: \"$name\""
   }
   if { $saveval ne "0"} {
    if { $saveval eq "-pruneby" || $saveval eq "-precision" || $saveval eq "-animate" } {
     set RADtoolflags "$RADtoolflags $saveval $nameo"
     unset singa($saveval)
    } else {
     set singa($saveval) $nameo
    }
    set saveval 0
    continue
   }
 
   if { [regexp -nocase {(^-bundle$)} $nameo] == 1 } {
    set saveflag 1
    continue
   } elseif { [regexp -nocase {(^-download$)} $nameo] == 1 } {
    set downloadflag 1
    continue
   } elseif { [regexp {(^-)} $nameo] == 1 } {
    if { [info exists singa($name)]} {
     # need to store the value.  For this, need to read next arg.
     set saveval $name
     continue
    }
    if { [regexp -nocase {(^-SSUonly$$)} $name] == 1 } {
     set asubunit 1
    }
    if { [regexp -nocase {(^-LSUonly$)} $name] == 1 } {
     set asubunit 2
    }
    set RADtoolflags "$RADtoolflags $nameo"
    set saveflag 0
    set downloadflag 0
    continue
   }
   if { $downloadflag == 1 } {
    set pdbdownload $nameo
    continue
   }
   if { $saveflag == 1 } {
    lappend pdblist $nameo
    continue
   }
   # if we made it to this point, then we are not reading pdb names, and it is not a standard flag.
   error "Argument format error: \"$nameo\""
  }
  if { [llength $pdblist] !=0 && [llength $pdbdownload]!=0 } {
   error "Can not use -bundle and -download flags together"
  } 
  if { [info exists bundleon] && [llength $pdblist]==0 } {
   error "A list of structure file names must be given after the -bundle flag"
  }
  if { [llength $pdbdownload]!=0 } {
   set pdblist [getpdb $pdbdownload]
  }
  foreach name [array names singa] { 
   if { $singa($name) eq "0"} {
    unset singa($name)  
   } 
  }
 
  return [list $pdblist $RADtoolflags [array get singa] $asubunit]
 }
 
 proc loadPDBs {pdblist} {
  puts "Input structure files (VMD molecule #, file name):"
  array set names {}
  foreach {name} $pdblist {
   set id [mol new "$name" waitfor -1]
   lappend IDs $id
   set names($id) $name
   puts "\t$id $name"
  }
  return [list [array get names] $IDs]
 }
 
 proc getCandidatesComplete {IDs} {
  set candidates {}
  puts "\nThe following chains are found to have at least 500\nnucleic acids and P atoms (VMD molecule #, chain ID, # P atoms)"
  foreach {molecule} $IDs {
   set sel [atomselect $molecule "name P and nucleic"]
   foreach {RNAchains} [lsort -unique [$sel get chain]] {
    set subsel [atomselect $molecule "chain \"$RNAchains\" and name P and nucleic"]
    set num [$subsel num]
    if {$num > 499} {
     lappend candidates [list $molecule $RNAchains $num]
     puts [list $molecule $RNAchains $num]
    }
    $subsel delete
   }
   $sel delete
  }
  # list them based on number of residues
  set candidates [lsort -integer -decreasing -index 2 $candidates]
  if { [llength $candidates] == 0 } {
   error "No chains were found that have more than 500 nucleic acids. 
 Sometimes this happens when RNA residues only contain P atoms.
 Nothing to do.  Will quit."
  }
  return $candidates
 }
 
 proc getCandidatesSingle {IDs SSULSU} {
  set candidates {}
  if { $SSULSU == 1 } {
   set name "SSU"
  } else {
   set name "LSU"
  }
  puts "\nThe following RNA molecules were found (VMD molecule #, chain ID, # P atoms)"
  foreach {molecule} $IDs {
   set sel [atomselect $molecule "name P and nucleic"]
   foreach {RNAchains} [lsort -unique [$sel get chain]] {
    set subsel [atomselect $molecule "chain \"$RNAchains\" and name P and nucleic"]
    set num [$subsel num]
    lappend candidates [list $molecule $RNAchains $num]
    puts [list $molecule $RNAchains $num]
    
    $subsel delete
   }
   $sel delete
  }
 
  puts "Will analyze chains with more than 500 nucleic acids:"
 
  set oldcandidates $candidates
  set candidates {}
  foreach {entry} $oldcandidates {
   if { [lindex $entry 2] >= 499} {
    puts [join $entry]
    lappend candidates $entry
   }
  }
 
  # list them based on number of residues
  set candidates [lsort -integer -decreasing -index 2 $candidates]
  return $candidates
 }
 
 proc runComplete {exec ofs nms candidates RADtoolflags} {
  array set names $nms
  array set outfiles $ofs
  set pairs {}
  array set paired {}
  set cands [llength $candidates]
  for {set i 0 } { $i < [expr $cands - 1] } {incr i} {
   if {! [info exists paired($i)]} {
    # not a paired subunit
    set maxconts 0
    set maxcontsID "NULL"
    set maxcontsID2 "NULL"
    set molID1 [lindex [lindex $candidates $i] 0]
    set chainID1 [lindex [lindex $candidates $i] 1]
    set sel1 [atomselect $molID1 "chain \"$chainID1\" and name P and nucleic"]
    for {set j [expr $i+1] } { $j < $cands } {incr j} {
     if {! [info exists paired($j)]} {
      set molID2 [lindex [lindex $candidates $j] 0]
      set chainID2 [lindex [lindex $candidates $j] 1]
      set sel2 [atomselect $molID2 "chain \"$chainID2\" and name P and nucleic"]
      set tin1 [lindex [$sel1 get index] 0]
      set tin1 [atomselect $molID1 "index $tin1"]
      set x1 [join [$tin1 get {x y z}]]
      $tin1 delete
      set tin2 [lindex [$sel2 get index] 0]
      set tin2 [atomselect $molID2 "index $tin2"]
      set x2 [join [$tin2 get {x y z}]]
      $tin2 delete
       
      if { [veclength [vecsub $x1 $x2] ] < 0.1  } {
       puts stderr "\n\n\nWarning: The first atoms in chain $chainID1 in file [molinfo $molID1 get filename] and chain $chainID2 in file [molinfo $molID2 get  filename] are less than 0.1 Angstroms apart. This is likely due to a mistake, such as giving two files that have the same structures. Most likely, RADtool will fail to complete with these files.\n\n"
       puts stderr "Type \"y\" if you want to continue, or else RADtool will exit."
       set resp [gets stdin]
       if { [regexp -nocase {(^[y]$)} $resp ] != 1 } {
        puts "Ok, will quit immediately..."
        $sel1 delete
        $sel2 delete
        return
       }
      }
      set contn [llength [lindex [ measure contacts 10 $sel1 $sel2 ] 0]]
      if {$contn > $maxconts} {
       set maxconts $contn
       set maxcontsID [list $names($molID2) $chainID2]
       set maxcontsID2 $j
       set maxmolID2 $molID2
      }
      $sel2 delete
     }
    }
    $sel1 delete
    if { $maxcontsID ne "NULL" } {
     puts "\nIdentified pair # $i:
 	LSU: $molID1, chain $chainID1
 	SSU: $maxmolID2, chain [lindex $maxcontsID 1]"
 
     lappend pairs [list $names($molID1) $chainID1] $maxcontsID
     set paired($maxcontsID2) 1
    }
   }
  }
  set i 0
  if {[llength $pairs] < 2 } {
   error "Was unable to automatically identify at least one LSU-SSU pair.
 This can happen for several reasons, such as:
      - the structure only contains a single subunit
      - you did not provide all structure files as input
      - This is an atypical ribosome that doesn't
          contain single large rRNAs in the LSU and SSU"
  }
 
  foreach {LSU SSU} $pairs {
   puts "\nBundle module will now call $exec to analyze LSU-SSU pair # $i"
   puts -nonewline "    Output files: "
   set IO {}
   foreach {name} [array names outfiles] {
    lappend IO $name
    if {[llength $pairs] > 2 } {
     puts -nonewline "$outfiles($name)\_$i "
     lappend IO "$outfiles($name)\_$i"
    } else {
     puts -nonewline "$outfiles($name) "
     lappend IO "$outfiles($name)"
    }
   }
   set first [expr  [molinfo num] ]

   set val 0 
   if {[catch {
    set val [$exec -l [lindex $LSU 0] -lc [lindex $LSU 1]  -s [lindex $SSU 0] -sc [lindex $SSU 1]  [join $RADtoolflags] [join $IO]] 
   } errorMessage] != 0 } {
    # this means stamp failed.  So, let's hide anything we just loaded.  We'll leave it in the 
    # session, in case someone wants to see what happened.
    set last [expr  [molinfo num] ]
    puts stderr "Attempt to analyze a pairs of chains failed. Will hide molecules loaded when trying (mol IDs [molinfo index $first] to [molinfo index [expr $last-1]]), though you may want to inspect them, in order to troubleshoot.\n\nError message returned by RADtool:\n$errorMessage"
    for {set j $first } { $j < $last } { incr j } {
     mol off [molinfo index $j] 
    }
   } else {
    set last [expr  [molinfo num] ]
    puts "\nAnalysis of the LSU-SSU pair completed. In the VMD main window, mol IDs [molinfo index $first] to [molinfo index [expr $last-1]] are associated with this assembly."
   }


   if { [info exists ::RADTOOL::radenv(ROTTEST)] } {
    return $val
   }
   incr i 
  }
 }
 
 proc runSingleSU {exec ofs nms candidates RADtoolflags asubunit} {
  array set names $nms
  array set outfiles $ofs
  if { $asubunit == 1 } {
   set flag1 "-s"
   set flag2 "-sc"
   set subunit "SSU"
  }
  if { $asubunit == 2 } {
   set flag1 "-l"
   set flag2 "-lc"
   set subunit "LSU"
  }
  set i 0
  foreach {SSU} $candidates {
   puts "\nWill now call $exec to analyze candidate $subunit # $i"
   puts -nonewline "    Output files:"
   set IO {}
   foreach {name} [array names outfiles] {
    lappend IO  $name
    if {[llength $candidates] > 1 } {
     puts -nonewline "$outfiles($name)\_$i "
     lappend IO "$outfiles($name)\_$i"
    } else {
     puts -nonewline "$outfiles($name) "
     lappend IO "$outfiles($name)"
    }
   }
   # note: we don't need to add -SSUonly, since it is carried through IO
   set first [molinfo num]
   set code 0 
   if {[catch {
    set code [$exec $flag1 $names([lindex $SSU 0]) $flag2 [lindex $SSU 1]  $RADtoolflags [join $IO]] 
   } errorMessage] != 0 } {
    # this means stamp failed.  So, let's hide anything we just loaded.  We'll leave it in the 
    # session, in case someone wants to see what happened.
    set last [expr  [molinfo num] ]
    puts stderr "Attempt to analyze a subunit failed. Will hide molecules loaded when trying (mol IDs [molinfo index $first] to [molinfo index [expr $last-1]]), though you may want to inspect them, in order to troubleshoot.\n\nError message returned by RADtool:$errorMessage\n"
    for {set j $first } { $j < $last } { incr j } {
     mol off [molinfo index $j] 
    }
   } else {
    set last [expr  [molinfo num] ]
    puts "\nAnalysis of the subunit completed. In the VMD main window, mol IDs [molinfo index $first] to [molinfo index [expr $last-1]] are associated with this subunit."
   }
   incr i 
  }
 }

}

