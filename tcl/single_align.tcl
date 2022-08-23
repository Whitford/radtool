###########################ORIGINAL COPYRIGHT#####################################
# University of Illinois Open Source License
# Copyright 2004-2007 Luthey-Schulten Group, 
# All rights reserved.
#
# $Id: multiseq.tcl,v 1.36 2015/05/20 22:19:08 kvandivo Exp $
# 
# Developed by: Luthey-Schulten Group
#      University of Illinois at Urbana-Champaign
#      http://www.scs.illinois.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
#
# Author(s): Elijah Roberts, John Eargle, Dan Wright, Michael Bach

# Versions
#   2.1 , revised to 3.0 - Kirby tweaks
#
####################END OF ORIGINAL COPYRIGHT#####################################

####################UPDATED/EXTENDED COPYRIGHT####################################
# This is a modified version of the method menu_stamp, which is part of multiseq. 
# It also contains an identical copy of concatenateRegionsIntoSequences
# While it could be possible to directly use these routines, as provided in 
# existing packages, duplicating the routines avoids possible variable collisions. 
# Modifications were introduced by the Whitford group, at Northeastern University. 
# Above is the original copyright information for multiseq. In accordance with the 
# above copyright, all terms apply to this derived version (single_align.tcl)   
# 
# The following addendum applies - Neither the names of the Whitford Group, 
# Northeastern University, nor the names of its contributors may be used to 
# endorse or promote products derived from this software without specific prior 
# written permission.
##################################################################################

namespace eval RTALIGN {
 
 
 # this performs a single sequence alignment with ClustalW
 proc single_seqalign {name molID1 chain1 residues1 molID2 chain2 residues2} {
 
  package require clustalw 1.1
 
  set ::ClustalW::tempDir $name
  if {[file exists $::ClustalW::tempDir]} {
   if {[catch {
    file delete -force $::ClustalW::tempDir
   } errorMessage] != 0} {
    puts "Temp directory $name already exists, and script was unable to remove it.  This is typically harmless, though you may end up with extra directories after using the script. Message returned when trying to remove directory:\n$errorMessage"
   }
  }
  if { ! [file exists $::ClustalW::tempDir]} {
   file mkdir $::ClustalW::tempDir
  }
 
  # this is a modified version of menu_stamp, provided by VMD.
  # major changes include:
  # only support for alignment of 2 regions
  # must give molIDs, chainIDs, and residue numbers (i.e. residue, not resid) of both regions
 
  # Set the stamp parameters: setting to default
  
  # Extract the appropriate regions.
  set preservePrefixSuffixAlignment 0
  
  #::SeqData::VMD::updateVMDSequences
  # first get the seqIDs from the molIDs and chains
  set seqID1 [addSequence $molID1 $chain1 ]
  set seqID2 [addSequence $molID2 $chain2 ]
  
  # if multiple sequences are associated with the selection, use the first.
  if {[llength $seqID1] > 1} {
   set seqID1 [lindex $seqID1 0]
  }
  if {[llength $seqID2] > 1} {
   set seqID2 [lindex $seqID2 0]
  }
  
  # want this 
  set seltoalign [list $seqID1 $seqID2  ]
  set alignedSequences [::ClustalW::alignSequences $seltoalign]
  if {[file exists $::ClustalW::tempDir]} {
   if {[catch {
    file delete -force $::ClustalW::tempDir
   } errorMessage] != 0} {
    puts "Temp directory $name already exists, and script was unable to remove it.  This is typically harmless, though you may end up with extra directories after using the script. Message returned when trying to remove directory:\n$errorMessage"
   }
  }
 
  return $alignedSequences
 }
 
 
 # -----------------------------------------------------------------------
 # This method performs a single structural alignment with stamp.
 # we could have probably directly use the stamp.tcl routines, but this works.
 proc single_stamp {name molID1 chain1 residues1 molID2 chain2 residues2} {
  package require stamp 1.2
  ::STAMP::setTempFileOptions $name $name
 
  if {[file exists $::STAMP::tempDir]} {
   if {[catch {
    file delete -force $::STAMP::tempDir
   } errorMessage] != 0} {
    puts "Temp directory $name already exists, and script was unable to remove it.  This is typically harmless, though you may end up with extra directories after using the script. Message returned when trying to remove directory:\n$errorMessage"
   }
  }
  if { ! [file exists $::STAMP::tempDir]} {
   file mkdir $::STAMP::tempDir
  }
 
  set sel1 [atomselect $molID1 "chain \"$chain1\" and residue $residues1 and name P"]
  set sel2 [atomselect $molID2 "chain \"$chain2\" and residue $residues2 and name P"]
  if { [$sel1 num] > [$sel2 num] } {
   # if more atoms in ID1, then switch order
   set tm $molID1
   set tc $chain1
   set tr $residues1
   set molID1 $molID2
   set chain1 $chain2
   set residues1 $residues2
   set molID2 $tm
   set chain2 $tc
   set residues2 $tr
   set rev 1
  } else {
   set rev 0
  }
  $sel1 delete
  $sel2 delete
  # this is a modified version of menu_stamp, provided by VMD.
  # major changes include:
  # only support for alignment of 2 regions
  # must give molIDs, chainIDs, and residue numbers (i.e. residue, not resid) of both regions
 
  # Set the stamp parameters: setting to default
  set npass 2
  set scanscore 6
  set scanslide 10
  # set to 1 to turn on
  set slowscan 0
  set scan 1
  
  # Extract the appropriate regions.
  set preservePrefixSuffixAlignment 0
  set extractedRegions {}
  
  #::SeqData::VMD::updateVMDSequences
  # first get the seqIDs from the molIDs and chains
  set seqID1 [addSequence $molID1 $chain1 ]
  set seqID2 [addSequence $molID2 $chain2 ]
 
  # get the first residue in each chain.
  set selt [atomselect $molID1 "chain \"$chain1\""]
  set firstres1 [lindex [$selt get residue] 0]
  $selt delete 
  set selt [atomselect $molID2 "chain \"$chain2\""]
  set firstres2 [lindex [$selt get residue] 0] 
  $selt delete 
  # then subtract the first residue from the residue numbers.  This is needed when constructing a list to replace the output of first residue in the chain
  set renum_residues1 {}
  foreach res $residues1 {
   set res [expr $res - $firstres1]
   lappend renum_residues1 $res	
  }
  set renum_residues2 {}
  foreach res $residues2 {
   set res [expr $res - $firstres2]
   lappend renum_residues2 $res	
  }
 
  set seltoalign [list $seqID1 $renum_residues1 $seqID2 $renum_residues2 ]
  # Extract the selected regions.
  set extractedRegions [::SeqData::extractRegionsFromSequences $seltoalign]
  if {$extractedRegions == {}} {
   error "The selection within each sequence must be contiguous."
  }
  # the rest of this is just the original menu_stamp method, except for the returned information.
  # Mark that we should preserve the rest of the alignment.
  set preservePrefixSuffixAlignment 1
  # Exract the region info.
  set originalSequenceIDs [lindex $extractedRegions 0]
  set regionSequenceIDs [lindex $extractedRegions 1]
  set prefixes [lindex $extractedRegions 2]
  set prefixEndPositions [lindex $extractedRegions 3]
  set suffixes [lindex $extractedRegions 4]
  set suffixStartPositions [lindex $extractedRegions 5]
  
  # Align the structures.
  if {[catch {
   set alignedSequenceIDs [::STAMP::alignStructures $regionSequenceIDs $scan $scanslide $scanscore $slowscan $npass]
  } errorMessage] != 0} {
   error "The STAMP alignment failed. STAMP errors are sometimes cryptic.  Please contact us <p.whitford@northeastern.edu> if you have questions. The following error message was returned by STAMP:\n$errorMessage"
  }
  # Put the sequences back together again.
  set originalSequenceIDs [concatenateRegionsIntoSequences $originalSequenceIDs $alignedSequenceIDs $prefixes $prefixEndPositions $suffixes $suffixStartPositions $preservePrefixSuffixAlignment]
  if {$originalSequenceIDs == {}} {
   error "Concatenation of the alignment data failed. The data may be in an inconsistent state."
  }
 
  if {[file exists $::STAMP::tempDir]} {
   if {[catch {
    foreach file [glob -nocomplain $::STAMP::tempDir/$::STAMP::filePrefix.*] {
     file delete -force $file
    }
    file delete -force $::STAMP::tempDir
   } errorMessage] != 0} {
    puts "Unable to remove temp directory $name.  This is typically harmless, though you may end up with extra directories after using the script. Message returned when trying to remove directory:\n$errorMessage"
   }
  }
 
  if { $rev == 1 } { 
   return [list $seqID2 $seqID1]
  } else {
   return [list $seqID1 $seqID2]
  }
 }
 
 # this is a modified version of updateVMDsequences.  It is used to load a single chain, rather than all chain sequences.  This was changed, so that we don't need to load every chain in a ribosome. Rather, we will just load the chain of interest. 
 
 proc addSequence { molID chainID } {
  # Go through each molecule in VMD.
  set addedSequenceIDs {}
 
  set chain $chainID
  set segname "" 
  # Get the atoms of the molecule.
  set atoms [atomselect $molID "chain \"$chain\" and not water and not ions"]
 
  # Make sure this parts is not a bunch of MD waters or lipids.
  set numatoms [$atoms num]
  if {$numatoms == 0} {
   error "Attempting to add a sequence with 0 atoms. MolID: $molID; Chain: $chain"
  }
  # Get the residue information.
  set atomPropsList [$atoms get {resid insertion \
                              resname residue altloc}]
  if {$numatoms == [llength $atomPropsList]} {
   # Go through each residue.
   set uniqueResIDs {}
   set uniqueResidues {}
   set residues {}
   set sequence {}
   foreach atomProps $atomPropsList {
 
    # Get the residue id.
    set resID [lindex $atomProps 0]
    set insertion [lindex $atomProps 1]
    set resName [lindex $atomProps 2]
    set residue [lindex $atomProps 3]
    set altloc [lindex $atomProps 4]
    # Make sure we process each unique resid only once.
    if {[lsearch $uniqueResIDs "$resID$insertion"]==-1} {
     lappend uniqueResIDs "$resID$insertion"
 
     # See if this is just a normal residue.
     if {$insertion == " " && ($altloc == "" || $altloc ==1) } {
      lappend residues $resID
      lappend sequence [::SeqData::VMD::lookupCode $resName]
     } else {
      # We must have either an insertion or 
      # an alternate location.
 
      # save the residue as a list.
      lappend residues [list $resID \
                             $insertion $altloc]
      lappend sequence [::SeqData::VMD::lookupCode $resName]
 
      # See if we should update any previous residues
      set updateIndex [lsearch $residues $resID]
      if {$updateIndex != -1} {
       set residues [lreplace $residues \
              $updateIndex $updateIndex \
                  [list $resID " " ""]]
      }
     }
    }          
 
    # Store the residue to check for VMD preload problem.
    if {[lsearch $uniqueResidues $residue] == -1} {
     lappend uniqueResidues $residue
    }                        
   }
 
   # If we only have one residue, make sure we only 
   # have one resid to avoid VMD preload problem.
   if {[llength $uniqueResidues] > 1 || \
                        [llength $uniqueResIDs] == 1} {
 
    # If there is more than one parts for this 
    # molecule, use the parts as the name suffix.
    set nameSuffix ""
 
    # Add the sequence to the sequence store.
    set sequenceID [::SeqData::VMD::addVMDSequence $molID \"$chain\" \
                $segname $sequence \
                "[::SeqData::VMD::getMoleculeName $molID]$nameSuffix" \
                [lindex $residues 0] $residues]
    if {$sequenceID != -1} {
     # Add sequence to list of seq ids we have added.
     lappend addedSequenceIDs $sequenceID
    }
 
   } else {
      break
   }
  }
  # Remove the atomselect.
  $atoms delete
  # Return the lists.
  if {[llength $addedSequenceIDs]>1} {
   error "When adding molID $molID, chain $chainID, multiple sequences returned"
  }
  return $addedSequenceIDs
 } 
 
 # This is simply the same proc as described in multiseq.tcl
 # it is copied here since it is the only thing we need from multiseq
 proc concatenateRegionsIntoSequences {originalSequenceIDs \
                              regionSequenceIDs prefixes prefixEndPositions \
                              suffixes suffixStartPositions \
                              preservePrefixSuffixAlignment} {
   
       set osiLength [llength $originalSequenceIDs]
       # Make sure all of the lists are the same length.
       if {[llength $regionSequenceIDs] != $osiLength ||\
             [llength $prefixes] != $osiLength || \
             [llength $prefixEndPositions] != $osiLength || \
             [llength $suffixes] != $osiLength || \
             [llength $suffixStartPositions] != $osiLength } {
          return {}
       }
 
       # Figure out the positioning.
       set maxPrefixEndPosition -1
       set minSuffixStartPosition -1
       for {set i 0} {$i < [llength $originalSequenceIDs]} {incr i} {
          set prefixEndPosition [lindex $prefixEndPositions $i]
          set suffix [lindex $suffixes $i]
          set suffixStartPosition [lindex $suffixStartPositions $i]
          if {$maxPrefixEndPosition == -1 || $prefixEndPosition > $maxPrefixEndPosition} {
             set maxPrefixEndPosition $prefixEndPosition
          }
          if {$suffix != {} && ($minSuffixStartPosition == -1 || $suffixStartPosition < $minSuffixStartPosition)} {
             set minSuffixStartPosition $suffixStartPosition
          }
       }
 
       # Copy the region sequence data into the original sequence ids.
       for {set i 0} {$i < [llength $originalSequenceIDs]} {incr i} {
 
          # Get the data for this sequence.
          set originalSequenceID [lindex $originalSequenceIDs $i]
          set regionSequenceID [lindex $regionSequenceIDs $i]
          set prefix [lindex $prefixes $i]
          set suffix [lindex $suffixes $i]
          set prefixEndPosition [lindex $prefixEndPositions $i]
          set suffixStartPosition [lindex $suffixStartPositions $i]
 
          # See if we are preserving the alignment of the prefixes and the suffixes.
          if {$preservePrefixSuffixAlignment} {
 
             # Get any gaps needed to keep the prefix in alignment.
             set prefixGaps [::SeqData::getGaps [expr $maxPrefixEndPosition-$prefixEndPosition]]
 
             # Get any gaps needed keep the suffix in alignment.
             set suffixGaps {}
             if {$suffix != {}} {
                set suffixGaps [::SeqData::getGaps [expr $suffixStartPosition-$minSuffixStartPosition]]
             }
 
             # Concatenate the sequence.
             ::SeqData::setSeq $originalSequenceID [concat $prefix $prefixGaps [::SeqData::getSeq $regionSequenceID] $suffixGaps $suffix]
 
          # Otherwise, ignore the alignment of the prefixes and suffixes.
          } else {
 
             # Get any gaps needed to gap out the prefix.
             set prefixGaps [::SeqData::getGaps [expr $maxPrefixEndPosition-$prefixEndPosition]]
 
             # Concatenate the sequence.
             ::SeqData::setSeq $originalSequenceID [concat $prefixGaps $prefix [::SeqData::getSeq $regionSequenceID] $suffix]
          }
       }
 
       return $originalSequenceIDs
    } ; # end of concatenateRegionsIntoSequences 
 
}  
