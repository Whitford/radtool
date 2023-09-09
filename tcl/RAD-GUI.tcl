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


# Declare global variables for this package.

package provide radgui 1.0

namespace eval RADGUI:: {

 # Export the package functions.
 namespace export RADGUI
 # main window 
 variable w
 # animation window
 variable wan
# testing a dialogue window
 variable RCSBname 
 variable LSUSSU 
 variable doSTAMP 
 variable doFINDHEAD2
 variable doDUMP 
 variable pruneON 
 variable pruneby 
 variable saveDir {} 
 variable savePrefix {} 
 variable animateON 
 variable LSUfile {}
 variable SSUfile {}
 variable LSUchain {}
 variable SSUchain {}
 variable stepnum {}
 variable startmol {}
 variable endmol {}
 variable availmods {}
 variable MODN {}
 variable VSEQ {}
 variable BUNDLEfile {}
 variable BUNDLEfilestr ""
 
 proc RAD-GUI {} {
  variable w
  variable RCSBname 
  variable LSUSSU 
  variable doSTAMP
  variable doFINDHEAD2
  variable doDUMP
  variable pruneON
  variable pruneby
  variable saveDir
  variable savePrefix
  variable animateON
  variable LSUfile 
  variable SSUfile 
  variable LSUchain 
  variable SSUchain 
  variable MODN {0}
  variable VSEQ {tr_r_t}
  variable BUNDLEfile {}
  variable BUNDLEfilestr ""

  set RCSBname {}
  set LSUSSU {all}
  set doSTAMP {1}
  set doFINDHEAD2 {0}
  set doDUMP {0}
  set pruneON {1}
  set pruneby {2}
  set animateON {none}
  
  set LSUfile {} 
  set SSUfile {}
  set LSUchain {}
  set SSUchain {}

  set saveDir {~/}
  set savePrefix {}

  if { [lsearch [font names] TkDefaultFontBold] != -1 } { 
   font delete TkDefaultFontBold
  }
  if { [lsearch [font names] TkFontBold2] != -1 } { 
   font delete TkFontBold2
  }
  set fontOpts {}
  foreach {flag val} [font configure TkDefaultFont] {lappend fontOpts $val}
  lassign $fontOpts font_f font_s font_w font_sl font_u font_o
  font create TkFontBold2 -family $font_f -size $font_s -weight bold -slant $font_sl -underline $font_u -overstrike $font_o 
  set font_s [expr {$font_s + 1}] 
  font create TkDefaultFontBold -family $font_f -size $font_s -weight bold -slant $font_sl -underline $font_u -overstrike $font_o 
  unset fontOpts font_f font_s font_w font_sl font_u font_o
 
  # create the dialogue
  # Find a name for the dialog
  set unique 0
  set prefix ".radtool"
  set childList [winfo children .]
  if {[lsearch $childList $prefix$unique] != -1} {
   wm deiconify $w
   return
  }
  set w [toplevel $prefix$unique]
  wm title $w "RADtool v$::RADTOOL::radversion"
  raise $w
  # Create the components.
  frame $w.all1
  frame $w.all1.structure 
  frame $w.all1.structure.a -relief sunken -borderwidth 1
  frame $w.all1.structure.a.b -relief raised -borderwidth 1
  label $w.all1.structure.a.b.head -text {Structure Information} -anchor center -font TkDefaultFontBold
  label $w.all1.structure.a.b.pdblabel1 -text "Option 1) Download and automatically identify the subunits" -anchor w -font TkFontBold2
  label $w.all1.structure.a.b.pdblabel2 -text {RCSB ID:} -anchor w
  entry $w.all1.structure.a.b.pdbname -width 4 -textvariable RADGUI::RCSBname

  label $w.all1.structure.a.b.head3 -text "Option 2) Provide files and detect LSU/SSU rRNA" -anchor w -font TkFontBold2 

  label $w.all1.structure.a.b.bundlelabel -text {Structure files:}
  label $w.all1.structure.a.b.bundleentry -width 20 -textvariable RADGUI::BUNDLEfilestr

  button $w.all1.structure.a.b.bundlebut -text {Add files} \
    -command {
      set tempfile [tk_getOpenFile -multiple 1 -parent .]
      if {![string equal $tempfile ""]} { 
       lappend RADGUI::BUNDLEfile $tempfile 
        set RADGUI::BUNDLEfilestr ""
        foreach filen [join $::RADGUI::BUNDLEfile] {
         set RADGUI::BUNDLEfilestr "$::RADGUI::BUNDLEfilestr[file tail $filen]\n"
        }
        set RADGUI::BUNDLEfilestr [string trim $::RADGUI::BUNDLEfilestr]
      }
    }
  button $w.all1.structure.a.b.bundlecbut -text {Clear} \
    -command {
     set RADGUI::BUNDLEfilestr ""
     set RADGUI::BUNDLEfile {} 
    }


  label $w.all1.structure.a.b.head2 -text "Option 3) Analyze specific files and chains" -anchor w -font TkFontBold2 

  label $w.all1.structure.a.b.lsulabel -text {File containing LSU:}
  entry $w.all1.structure.a.b.lsuentry -width 20 -textvariable RADGUI::LSUfile

  button $w.all1.structure.a.b.lsubut -text {Browse} \
    -command {
      set tempfile [tk_getOpenFile -parent .]
      if {![string equal $tempfile ""]} { set RADGUI::LSUfile $tempfile }
    }

  label $w.all1.structure.a.b.lsuclabel -text {LSU chain ID:}
  entry $w.all1.structure.a.b.lsucentry -width 4 -textvariable RADGUI::LSUchain

  label $w.all1.structure.a.b.ssulabel -text {File containing SSU:}
  entry $w.all1.structure.a.b.ssuentry -width 20 -textvariable RADGUI::SSUfile

  button $w.all1.structure.a.b.ssubut -text {Browse} \
    -command {
      set tempfile [tk_getOpenFile -parent .]
      if {![string equal $tempfile ""]} { set RADGUI::SSUfile $tempfile }
    }

  label $w.all1.structure.a.b.ssuclabel -text {SSU chain ID:}
  entry $w.all1.structure.a.b.ssucentry -width 4 -textvariable RADGUI::SSUchain

  frame $w.all2
  frame $w.all2.analysis
  label $w.all2.analysis.head -text {Analysis and Animation Options} -anchor center -font TkDefaultFontBold
  label $w.all2.analysis.choice -text {Subunits to analyze:}
  frame $w.all2.analysis.g1
   frame $w.all2.analysis.g1.a -relief sunken -borderwidth 1
   frame $w.all2.analysis.g1.a.b -relief raised -borderwidth 1
   radiobutton $w.all2.analysis.g1.a.b.lsussu -text {LSU-SSU pairs} -variable RADGUI::LSUSSU -value "all"
   radiobutton $w.all2.analysis.g1.a.b.ssuonly -text {SSUs only} -variable RADGUI::LSUSSU -value "-SSUonly"
   radiobutton $w.all2.analysis.g1.a.b.lsuonly -text {LSUs only} -variable RADGUI::LSUSSU -value "-LSUonly"

  label $w.all2.analysis.stampLabel -text {Perform STAMP alignment:}
  checkbutton $w.all2.analysis.stamp -variable RADGUI::doSTAMP -onvalue 1 -offvalue 0
  
  label $w.all2.analysis.pruneLabel -text {Prune cores:}
  checkbutton $w.all2.analysis.prune -variable RADGUI::pruneON -onvalue 1 -offvalue 0

  label $w.all2.analysis.fhLabel -text {Alt. find head method:}
  checkbutton $w.all2.analysis.fh -variable RADGUI::doFINDHEAD2 -onvalue 1 -offvalue 0

  label $w.all2.analysis.prunebyLabel -text {Prune threshold (Angstroms):} -anchor w
  entry $w.all2.analysis.pruneby -width 4 -textvariable RADGUI::pruneby  
 
  label $w.all2.analysis.choice2 -justify left -text {Animation sequence:
Will show the reference/
classical model (4V9D),
and move it to the model(s)
being analyzed. This can
be memory/RAM intensive.} 
  frame $w.all2.analysis.g2
   frame $w.all2.analysis.g2.a -relief sunken -borderwidth 1
   frame $w.all2.analysis.g2.a.b -relief raised -borderwidth 1
   radiobutton $w.all2.analysis.g2.a.b.off -text   { none  } -variable RADGUI::animateON -value "none"
   radiobutton $w.all2.analysis.g2.a.b.trrot -text { translate, rotate, tilt } -variable RADGUI::animateON -value "tr_r_t"
   radiobutton $w.all2.analysis.g2.a.b.trtro -text { translate, tilt, rotate } -variable RADGUI::animateON -value "tr_t_r"
   radiobutton $w.all2.analysis.g2.a.b.rottr -text { rotate, tilt, translate } -variable RADGUI::animateON -value "r_t_tr"
   radiobutton $w.all2.analysis.g2.a.b.trotr -text { tilt, rotate, translate } -variable RADGUI::animateON -value "t_r_tr"

  frame $w.all3
  frame $w.all3.output
  frame $w.all3.output.a -relief sunken -borderwidth 1
  frame $w.all3.output.a.b -relief raised -borderwidth 1
  label $w.all3.output.a.b.head -text {Output Options} -anchor center -font TkDefaultFontBold
  label $w.all3.output.a.b.saveInfoLab -justify left -text {If you want to save output to files, then specify a prefix and directory.
Notes: If more than one ribosome is analyzed with options 1 or 2, then
an integer will be inserted before the suffix of each output file name.
Downloaded files will be saved to the same directory as output.}
  label $w.all3.output.a.b.saveDirLab -text {Select output directory:}
  entry $w.all3.output.a.b.saveDirEntry -width 20 -textvariable RADGUI::saveDir

  button $w.all3.output.a.b.saveDirBut -text {Browse} \
   -command {
   RADGUI::chooseSaveDir
  }
   # this was modeled after various examples in VMD tcl scripts 
   proc chooseSaveDir {} {
    set gooddir 0
    while { $gooddir == 0 } {
     set dir [ tk_chooseDirectory -mustexist true -title "Choose Save Directory" -initialdir RADGUI::saveDir  ]
     if { $dir == "" } {
      return
     }
     if { [file writable "$dir"] } {
      set gooddir 1
      set RADGUI::saveDir $dir
     } else {
      tk_messageBox -type ok -message "Error: $dir is not writable"
     }
    }
   }
     
  label $w.all3.output.a.b.savePrefixLab   -text {Output file prefix:}
  entry $w.all3.output.a.b.savePrefixEntry -width 20 -textvariable RADGUI::savePrefix

  label $w.all3.output.a.b.dumpLabel -text {Print a lot of stuff:}
  checkbutton $w.all3.output.a.b.dump -variable RADGUI::doDUMP -onvalue 1 -offvalue 0

  frame $w.all3.opts
  button $w.all3.opts.http -text "RADtool webpage" -command {vmd_open_url http://www.radtool.org}
  button $w.all3.opts.details -text "More Information" -command {::RADGUI::moreinfo}
  button $w.all3.opts.animate -text "More Animation Options" -command {::RADGUI::animate2}

  frame $w.all3.bottom
   frame $w.all3.bottom.buttons
    button $w.all3.bottom.buttons.accept   -text {ANALYZE STRUCTURE} -pady 2 -command RADGUI::but_run
    button $w.all3.bottom.buttons.clear    -text {RESET FORM} -pady 2 -command RADGUI::reset
    button $w.all3.bottom.buttons.cancel   -text {CLOSE} -pady 2 -command RADGUI::but_cancel
  frame $w.all3.contact
  label $w.all3.contact.info -justify left -text {RADtool is an external plugin developed by the Whitford Group at the
Center for Theoretical Biological Physics, Northeastern University.
For support, contact us and/or consult the RADtool webpage.}
 
  # Layout the components.
  grid $w.all1                               -row 1 -column 1 -sticky n
  pack $w.all1.structure                     -fill both -expand true -side top -padx 5 -pady 5
  pack $w.all1.structure.a                   -fill both -expand true -side left
  pack $w.all1.structure.a.b                 -fill both -expand true -side left
  grid $w.all1.structure.a.b.head            -row 1 -column 1 -columnspan 3 -sticky nw   
  grid $w.all1.structure.a.b.pdblabel1       -row 2 -column 1 -columnspan 3 -sticky w
  grid $w.all1.structure.a.b.pdblabel2       -row 3 -column 1 -columnspan 1 -sticky w
  grid $w.all1.structure.a.b.pdbname         -row 3 -column 2 -columnspan 1 -sticky w
  grid $w.all1.structure.a.b.head3           -row 4 -column 1 -columnspan 2 -sticky w 
  grid $w.all1.structure.a.b.bundlelabel     -row 5 -column 1 -sticky w
  grid $w.all1.structure.a.b.bundleentry     -row 5 -column 2 -sticky w
  grid $w.all1.structure.a.b.bundlebut       -row 5 -column 3 -sticky w
  grid $w.all1.structure.a.b.bundlecbut      -row 6 -column 3 -sticky w
  grid $w.all1.structure.a.b.head2           -row 7 -column 1 -columnspan 2 -sticky w 
  grid $w.all1.structure.a.b.lsulabel        -row 8 -column 1 -sticky w
  grid $w.all1.structure.a.b.lsuentry        -row 8 -column 2 -sticky w
  grid $w.all1.structure.a.b.lsubut          -row 8 -column 3 -sticky w
  grid $w.all1.structure.a.b.lsuclabel       -row 9 -column 1 -sticky w
  grid $w.all1.structure.a.b.lsucentry       -row 9 -column 2 -sticky w
  grid $w.all1.structure.a.b.ssulabel        -row 10 -column 1 -sticky w
  grid $w.all1.structure.a.b.ssuentry        -row 10 -column 2 -sticky w
  grid $w.all1.structure.a.b.ssubut          -row 10 -column 3 -sticky w
  grid $w.all1.structure.a.b.ssuclabel       -row 11 -column 1 -sticky w
  grid $w.all1.structure.a.b.ssucentry       -row 11 -column 2 -sticky w
 
 
  grid $w.all2                               -row 1 -column 2 -sticky n

  pack $w.all2.analysis                  -fill both -expand true -side top -padx 5 -pady 5
  grid $w.all2.analysis.head             -column 1 -row 1 -columnspan 3 -sticky nw 
  grid $w.all2.analysis.choice           -column 1 -row 2 -sticky nw -pady 3
  grid $w.all2.analysis.g1               -column 2 -row 2 -sticky nw -padx 5 -pady 3 -columnspan 2
  pack $w.all2.analysis.g1.a             -fill both -expand true -side left
  pack $w.all2.analysis.g1.a.b           -fill both -expand true -side left
  grid $w.all2.analysis.g1.a.b.lsussu    -column 1 -row 1 -sticky w
  grid $w.all2.analysis.g1.a.b.ssuonly   -column 1 -row 2 -sticky w
  grid $w.all2.analysis.g1.a.b.lsuonly   -column 1 -row 3 -sticky w
  grid $w.all2.analysis.stampLabel       -column 1 -row 3 -sticky nw -pady 3
  grid $w.all2.analysis.stamp            -column 2 -row 3 -sticky nw -padx 1 -pady 3 -columnspan 2
  grid $w.all2.analysis.fhLabel          -column 1 -row 4 -sticky nw -pady 3
  grid $w.all2.analysis.fh               -column 2 -row 4 -sticky nw -padx 1 -pady 3 -columnspan 2
  grid $w.all2.analysis.pruneLabel       -column 1 -row 5 -sticky nw -pady 3
  grid $w.all2.analysis.prune            -column 2 -row 5 -sticky nw -padx 1 -pady 3 -columnspan 2
  grid $w.all2.analysis.prunebyLabel     -column 1 -row 6 -sticky nw -pady 3
  grid $w.all2.analysis.pruneby          -column 2 -row 6 -sticky nw -padx 1 -pady 3 -columnspan 2

  grid $w.all2.analysis.choice2          -column 1 -row 7 -sticky nw -pady 3
  grid $w.all2.analysis.g2               -column 2 -row 7 -sticky nw -padx 5 -pady 3 -columnspan 2
  pack $w.all2.analysis.g2.a             -fill both -expand true -side left
  pack $w.all2.analysis.g2.a.b           -fill both -expand true -side left
  grid $w.all2.analysis.g2.a.b.off     -column 1 -row 1 -sticky w
  grid $w.all2.analysis.g2.a.b.trrot   -column 1 -row 2 -sticky w
  grid $w.all2.analysis.g2.a.b.trtro   -column 1 -row 3 -sticky w
  grid $w.all2.analysis.g2.a.b.rottr   -column 1 -row 4 -sticky w
  grid $w.all2.analysis.g2.a.b.trotr   -column 1 -row 5 -sticky w

  grid $w.all3                               -row 1 -column 3 -sticky n

  pack $w.all3.output                      -fill both -expand true -side top -padx 5 -pady 5
  pack $w.all3.output.a                    -fill both -expand true -side left
  pack $w.all3.output.a.b                  -fill both -expand true -side left
  grid $w.all3.output.a.b.head             -column 1 -row 1 -columnspan 3 -sticky nw 
  grid $w.all3.output.a.b.saveInfoLab      -row 2 -column 1 -columnspan 3 -sticky w
  grid $w.all3.output.a.b.savePrefixLab    -row 3 -column 1 -sticky w
  grid $w.all3.output.a.b.savePrefixEntry  -row 3 -column 2
  grid $w.all3.output.a.b.saveDirLab       -row 4 -column 1 -sticky w
  grid $w.all3.output.a.b.saveDirEntry     -row 4 -column 2
  grid $w.all3.output.a.b.saveDirBut       -row 4 -column 3
  grid $w.all3.output.a.b.dumpLabel            -column 1 -row 5 -sticky nw -pady 3
  grid $w.all3.output.a.b.dump                 -column 2 -row 5 -sticky nw -padx 1 -pady 3 -columnspan 2

  pack $w.all3.opts                    -fill x -side top
  pack $w.all3.opts.http               -side right  
  pack $w.all3.opts.details            -side left  
  pack $w.all3.opts.animate            -side left  
 
  pack $w.all3.bottom                  -fill x 
  pack $w.all3.bottom.buttons          
  pack $w.all3.bottom.buttons.accept   -side left -padx 5 -pady 5
  pack $w.all3.bottom.buttons.clear    -side left -pady 5
  pack $w.all3.bottom.buttons.cancel   -side left -pady 5
  pack $w.all3.contact                  -fill x -side bottom
  pack $w.all3.contact.info          -side bottom

  bind $w <Destroy> {RADGUI::but_cancel}
 
  if { [vmdinfo version] ne "1.9.3" } {
   if { [vmdinfo version] eq "1.9.4" } {
    tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "WARNING: YOU ARE USING VMD 1.9.4. THIS VERSION OF VMD WAS NOT AVAILABLE AT THE TIME RADTOOL WAS RELEASED. SO, THERE IS NO WAY TO KNOW IF THE PERFORMANCE WILL BE RELIABLE. CHECK AND SEE IF A NEW VERSION OF RADTOOL IS AVAILABLE, OR PERHAPS USE VMD 1.9.3."
    return
   } else {
tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "WARNING: YOU ARE USING VMD [vmdinfo version]!!!!
 RADtool HAS ONLY BEEN TESTED WITH VMD v1.9.3. PERFORMANCE WITH PRE-RELEASE VERSIONS OF VMD MAY BE LESS STABLE.
 
IF YOU RUN \"RADtool -test\" IN THE TkConsole, AND ALL TESTS ARE PASSED, THEN MOST FEATURES MAY BE WORKING PROPERLY. HOWEVER, SEE README FOR KNOWN ISSUES."
    return
   }
  }
 
 }

 proc moreinfo {} {
  variable w 
  menu tkcon on
  ::RADTOOL::splashmessage "stdout"
  ::RADTOOL::description
    return 

 }
 proc but_cancel {} {
  variable w
  # hide the dialog.
  wm withdraw $w
  reset
 }

 proc but_run {} {
  variable w
  variable RCSBname
  variable LSUSSU
  variable LSUonly
  variable SSUonly
  variable doSTAMP
  variable doFINDHEAD2
  variable doDUMP
  variable pruneON
  variable pruneby
  variable saveDir
  variable savePrefix
  variable LSUfile
  variable SSUfile
  variable LSUchain
  variable SSUchain
  variable animateON
  variable BUNDLEfile
  set origDir {}
  if { $RCSBname eq {} && $LSUfile eq {} && $SSUfile eq {}  && $BUNDLEfile eq {}} {
   tk_messageBox -icon error -type ok -title Message -parent $w \
   -message "Must give an RCSB ID or provide files for analysis!"
   return 
  }
  if { $RCSBname ne {} } {
   if {$LSUfile ne {} || $SSUfile ne {} || $LSUchain ne {} || $SSUchain ne {} || $BUNDLEfile ne {}} {
   tk_messageBox -icon error -type ok -title Message -parent $w \
   -message "Can not give RCSB entry along with specific structure files or chain IDs!"
   return
   } 
   if { [string length $RCSBname] != 4 } {
    tk_messageBox -icon error -type ok -title Message -parent $w \
    -message "Must give 4-character RCSB ID."
    return 
   } 
  } elseif {$BUNDLEfile ne {}} {
   if {$LSUfile ne {} || $SSUfile ne {} || $LSUchain ne {} || $SSUchain ne {}} {
   tk_messageBox -icon error -type ok -title Message -parent $w \
   -message "Can not give list of files (option 2), as well as specific LSU/SSU files or chain IDs (option 3)!"
   return
   } 
  } 
 
  set radmess {}
  set args {}
  if { $RCSBname ne {} } {
   # download the system
   lappend args "-download"
   lappend args  $RCSBname
   lappend radmess "Will try to download RCSB entry $RCSBname.\n"
   if { $saveDir eq {} } {
    tk_messageBox -icon error -type ok -title Message -parent $w \
    -message "When using download option, must specify directory for saving output."
    return 
   }
   if {[file exists $saveDir] } {
    if {![file isdirectory $saveDir]} {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "$saveDir does not appear to be a directory."
     return 
    }
    if {![file writable $saveDir]} {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Output directory ($saveDir) is not writable.\nPlease select a different output directory."
     return
    }
   } else {
    if { [catch {file mkdir $saveDir } errorMessage] != 0} {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Unable to make directory $saveDir. Perhaps you don\'t have the proper permissions to create this directory."
     return 
    }
   }
  } elseif {$BUNDLEfile ne {}} {
   lappend args "-bundle"
   lappend radmess "     Will read files and try to automatically detect LSUs and/or SSUs.\n      Will read:\n"
   foreach name [join $BUNDLEfile] {
    lappend radmess "      $name\n"
    lappend args $name
   } 
  } else {
   # must give files, since we are not downloading
   if { $LSUSSU eq "all" || $LSUSSU eq "-LSUonly" } {
    if { $LSUfile eq {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Must specify LSU file name."
     return 
    } else {
     lappend args "-l"
     lappend args $LSUfile
     lappend radmess "    Will read LSU file $LSUfile.\n"
    }
    if { $LSUchain eq {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Must specify LSU chain ID."
     return 
    } else {
     lappend args "-lc"
     lappend args $LSUchain
     lappend radmess "    Will define LSU chain ID as $LSUchain.\n"
    }
   }

   if { $LSUSSU eq "all" || $LSUSSU eq "-SSUonly" } {
    if { $SSUfile eq {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Must specify SSU file name."
     return 
    } else {
     lappend args "-s"
     lappend args $SSUfile
     lappend radmess "    Will read LSU file $SSUfile.\n"
    }
    if { $SSUchain eq {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Must specify SSU chain ID."
     return 
    } else {
     lappend args "-sc"
     lappend args $SSUchain
     lappend radmess "    Will define SSU chain ID as $SSUchain.\n"
    }
   }
  }

  if { $LSUSSU ne "all" } {
   lappend args $LSUSSU
    lappend radmess "    Will try to analyze LSU-SSU pairs.\n"
   if { $LSUSSU eq "-LSUonly" } {
    lappend radmess "    Will try to align LSUs, only.\n"
    if { $SSUfile ne {} || $SSUchain ne {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Can not give SSU information when analyzing LSUs only."
     return 
    }
   }
   if { $LSUSSU eq "-SSUonly" } {
    lappend radmess "    Will try to analyze SSUs, only.\n"
    if { $LSUfile ne {} || $LSUchain ne {} } {
     tk_messageBox -icon error -type ok -title Message -parent $w \
     -message "Can not give LSU information when analyzing SSUs only."
     return 
    }

   }
  }
  if { $doSTAMP eq 1 } {
   lappend args "-stamp"
   lappend radmess "    STAMP alignment will be performed.\n"
  }
  if { $doFINDHEAD2 eq 1 } {
   lappend args "-findhead2"
   lappend radmess "    Alternate method for finding head will be used.\n"
  }
  if { $doDUMP eq 1 } {
   lappend args "-dump"
   lappend radmess "    Will write lots of data (not usually needed).\n"
  }
 
  if { $animateON ne "none" } {
   lappend args "-animate"
   lappend args "$animateON"
   lappend radmess "    Turning on animation. NOTE: An animation will be created for each ribosome being analyzed.  So, if you have a bundle with multiple ribosomes, there will be multiple animations created.\n"
  }
  
  if { $pruneON eq 1 } {
   if { ! [ string is double $pruneby] && ! [ string is integer $pruneby] } {
    tk_messageBox -icon error -type ok -title Message -parent $w \
    -message "Pruning threshold value must be an integer or double.  Found \"$pruneby\""
    return 
   }
   lappend args "-pruneby"
   lappend args "$pruneby"
   lappend radmess "    Pruning will be performed, with a threshold of $pruneby Angstroms.\n"
  }
  if { $savePrefix ne {} } {
   if { $saveDir eq {} } {
    tk_messageBox -icon error -type ok -title Message -parent $w \
    -message "Save directory must be specified if prefix is given."
    return 
   }
   lappend radmess "    Will save files with prefix \"$savePrefix\"\n" 
   lappend radmess "    Will save files in directory \"$saveDir\"\n"
   lappend args "-o" 
   lappend args "$savePrefix.out" 
   lappend args "-e" 
   lappend args "$savePrefix.stderr" 
   lappend args "-cores_out" 
   lappend args "$savePrefix.cores" 
   if { $doSTAMP eq 1 } {
    lappend args "-align_out" 
    lappend args "$savePrefix.alignments" 
   }
  }
  puts "\n\nRAD-GUI STARTING\n[join $radmess]\n"
  if { $saveDir ne {} } {
   # we are just going to move to the save directory, so that all data will be dumped in the same place
   set origDir [pwd]
   cd $saveDir
  }
  set tkmess 0
  set err {}
  menu tkcon on
  if {[catch {
   set ::RADTOOL::radenv(GUION) 0 
   set err [::RADTOOL::RADtool $args]
   unset ::RADTOOL::radenv(GUION)  
  } errorMessage] != 0} {
   puts stderr "\nRADtool exited with the following error:\n\n$errorMessage"
   set tkmess 1
  }
  if {$err ne {}} {
   set tkmess 1
  }
  if {$tkmess == 1} {
   tk_messageBox -icon error -type ok -title Message -parent $w \
    -message "RADtool exited with an error. Check TkConsole and/or output files for details."
    return 
  }
  puts "\nRAD-GUI COMPLETED\n"
  if { $origDir ne {} } {
   # we are just going to move to the save directory, so that all data will be dumped in the same place
   cd $origDir 
  }

 }

 proc animate2 {} {
  variable wan
  variable stepnum {10}
  variable startmol {}
  variable endmol {}
  variable numIDs 
  set unique 0
  set prefix ".radtoolanimate"
  set childList [winfo children .]
  if {[lsearch $childList $prefix$unique] != -1} {
   #dialogue already open, so refresh

   but_refresh_an 
   return
  }

  # define this variable after the above lines, so that we only clear the variable when creating a new dialogue
  set wan [toplevel $prefix$unique]
  wm title $wan "RADtool Animations"
# make the list of options and add them to this label.
  set lb "#DDFFFE"
  set ly "#FDFFD5"
  set llg "#F8F8F8"
  set lg "#E7E7E7"

  makelist

  if { $numIDs == 0 } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "No LSU-SSU pairs or SSUs\n available for visualization.\nMust analyze some\nstructures, first."
   RADGUI::but_cancel_an
   return
  }

  frame $wan.opts 
  label $wan.opts.select -text {Animation Settings} -anchor center -font TkDefaultFontBold
  label $wan.opts.desc -text {This tool will generate an animation
between the orientations defined by
any two models that have been 
analyzed in the current session.
Note that, if a model or elements
of the model are deleted in the 
VMD window, it will not be shown here.
The Euler-Rodrigues angle between
models will also be calculated.
Select which models should be used
to define the starting and ending
orientations for the animation. 
Also, select which structure should
be displayed.} 
  label $wan.opts.startlabel -text {Starting Model #}
  label $wan.opts.endlabel -text {Ending Model #}
  label $wan.opts.numlabel -text "Frames per motion\ne.g. Rot, Tilt, Trans"
  entry $wan.opts.startmol -width 3 -textvariable ::RADGUI::startmol
  entry $wan.opts.endmol -width 3 -textvariable ::RADGUI::endmol
  entry $wan.opts.num -width 3 -textvariable ::RADGUI::stepnum
  label $wan.opts.choice -text {In animation, show:}
  frame $wan.opts.g1
   frame $wan.opts.g1.a -relief sunken -borderwidth 1
   frame $wan.opts.g1.a.b -relief raised -borderwidth 1
   radiobutton $wan.opts.g1.a.b.rbmod -text {Classical Model} -variable RADGUI::MODN -value 0
   radiobutton $wan.opts.g1.a.b.mod -text {Starting Model} -variable RADGUI::MODN -value 1
  label $wan.opts.vizseq -text {Animation sequence:}
  frame $wan.opts.g2
   frame $wan.opts.g2.a -relief sunken -borderwidth 1
   frame $wan.opts.g2.a.b -relief raised -borderwidth 1
   radiobutton $wan.opts.g2.a.b.sTrRT -text {trans, rot, tilt} -variable RADGUI::VSEQ -value tr_r_t
   radiobutton $wan.opts.g2.a.b.sTrTR -text {trans, tilt, rot} -variable RADGUI::VSEQ -value tr_t_r
   radiobutton $wan.opts.g2.a.b.sRTTr -text {rot, tilt, trans} -variable RADGUI::VSEQ -value r_t_tr
   radiobutton $wan.opts.g2.a.b.sTRTr -text {tilt, rot, trans} -variable RADGUI::VSEQ -value t_r_tr

  frame $wan.bottom
   frame $wan.bottom.buttons
    button $wan.bottom.buttons.accept   -text {Animate} -pady 2 -command RADGUI::but_run_an
    button $wan.bottom.buttons.refresh   -text {Refresh} -pady 2 -command RADGUI::but_refresh_an
    button $wan.bottom.buttons.cancel   -text {Close} -pady 2 -command RADGUI::but_cancel_an

  packlist  

  pack $wan.opts
  grid $wan.opts.select            -row 1 -column 1 -columnspan 3  -sticky news   
  grid $wan.opts.desc              -row 2 -column 1 -columnspan 3 -sticky news   
  grid $wan.opts.startlabel        -row 3 -column 1 -columnspan 1 -sticky w   
  grid $wan.opts.startmol          -row 3 -column 2 -columnspan 1 -sticky w   
  grid $wan.opts.endlabel          -row 4 -column 1 -columnspan 1 -sticky w   
  grid $wan.opts.endmol            -row 4 -column 2 -columnspan 1 -sticky w   
  grid $wan.opts.numlabel          -row 5 -column 1 -columnspan 1 -sticky w   
  grid $wan.opts.num               -row 5 -column 2 -columnspan 1 -sticky w   
  grid $wan.opts.choice           -column 1 -row 6 -sticky nw -pady 3
  grid $wan.opts.g1               -column 2 -row 6 -sticky nw -padx 5 -pady 3 -columnspan 2
  pack $wan.opts.g1.a             -fill both -expand true -side left
  pack $wan.opts.g1.a.b           -fill both -expand true -side left
  grid $wan.opts.g1.a.b.rbmod     -column 1 -row 1 -sticky w
  grid $wan.opts.g1.a.b.mod       -column 1 -row 2 -sticky w

  grid $wan.opts.vizseq           -column 1 -row 7 -sticky nw -pady 3
  grid $wan.opts.g2               -column 2 -row 7 -sticky nw -padx 5 -pady 3 -columnspan 2
  pack $wan.opts.g2.a             -fill both -expand true -side left
  pack $wan.opts.g2.a.b           -fill both -expand true -side left
  grid $wan.opts.g2.a.b.sTrRT       -column 1 -row 1 -sticky w
  grid $wan.opts.g2.a.b.sTrTR       -column 1 -row 2 -sticky w
  grid $wan.opts.g2.a.b.sRTTr       -column 1 -row 3 -sticky w
  grid $wan.opts.g2.a.b.sTRTr       -column 1 -row 4 -sticky w

  pack $wan.bottom                  -fill x -side bottom
  pack $wan.bottom.buttons          -side bottom
  pack $wan.bottom.buttons.accept   -side left -padx 5 -pady 5
  pack $wan.bottom.buttons.refresh  -side left -padx 5 -pady 5
  pack $wan.bottom.buttons.cancel   -side left -pady 5

  bind $wan <Destroy> {::RADGUI::but_cancel_an}
 }

 proc makelist {} {
  variable wan
  variable availmods {}
  variable numIDs {0}
  set lb "#DDFFFE"
  set ly "#FDFFD5"
  set llg "#F8F8F8"
  set lg "#E7E7E7"

  set RTHIST $::RADTOOL::radenv(RADtoolhistory)
  frame $wan.list 
  label $wan.list.label -text {Analyzed Models} -anchor center -font TkDefaultFontBold
  label $wan.list.label0 -text {Model #} -anchor center -font TkFontBold2 -bg $lg
  label $wan.list.label1 -text {SSU ID} -anchor center  -font TkFontBold2 -bg $llg
  label $wan.list.label2 -text {LSU ID} -anchor center  -font TkFontBold2 -bg $lg
  label $wan.list.summary -text "Model # is the index used by RADtool.\nSSU ID and LSU ID are the ID numbers\nlisted in the \"VMD Main\" window." -anchor center 

  set ii 0
  for {set i 0 } { $i < [llength $RTHIST]} { incr i} {
   set entry [lindex $RTHIST $i]
   set LSUID [lindex $entry end-1]
   set SSUID [lindex $entry end]
   set mod1 [lindex $entry end-4]
   set mod2 [lindex $entry end-3]
   if { (([lsearch [molinfo list] $LSUID] != -1 &&  [lsearch [molinfo list] $mod1] != -1  ) || $LSUID == -100 ) && [lsearch [molinfo list] $SSUID] != -1 && [lsearch [molinfo list] $mod2] != -1 } {
    if { [expr $ii % 2] == 0 } {
     set c1 $llg
     set c2 $lg
    } else {
     set c1 $lg
     set c2 $llg
    }
    label $wan.list.i$ii -text $i -bg $c1
    label $wan.list.s$ii -text $SSUID -bg $c2
    if {$LSUID == -100 } {
     set LSUID "N/A"
    }
    label $wan.list.l$ii -text $LSUID -bg $c1
    lappend availmods $i
    incr ii
   }
  }
  set numIDs $ii
 }

 proc packlist {} {
  variable wan
  variable numIDs
  pack $wan.list                     -fill both -expand true -side top -padx 5 -pady 5
  grid $wan.list.label            -row 1 -column 1 -columnspan 3 -sticky n   
  grid $wan.list.label0            -row 2 -column 1  -sticky news   
  grid $wan.list.label1            -row 2 -column 2  -sticky news   
  grid $wan.list.label2            -row 2 -column 3  -sticky news   
  for {set i 0 } { $i < $numIDs} { incr i} {
   grid $wan.list.i$i              -row [expr $i+3] -column 1 -sticky news 
   grid $wan.list.s$i              -row [expr $i+3] -column 2 -sticky news 
   grid $wan.list.l$i              -row [expr $i+3] -column 3 -sticky news 
  }
  grid $wan.list.summary            -row [expr $numIDs+3] -column 1 -columnspan 3  -sticky nw   
 }

 proc but_cancel_an {} {
  variable wan
  # Destroy the dialog.
  catch {
   bind $wan <Destroy> {}
   destroy $wan
  }
 }

 proc but_refresh_an {} {
  but_cancel_an
  animate2
 }

 proc but_run_an {} {
  variable wan
  variable stepnum
  variable startmol
  variable endmol
  variable availmods
  variable MODN
  variable VSEQ
  if { $stepnum eq {} } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Must specify number of frames per motion."
   return
  }
  if { $startmol eq {} } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Must specify starting model."
   return
  }
  if { $endmol eq {} } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Must specify ending model."
   return
  }
  if {[lsearch $availmods $startmol] == -1 } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Starting model \"$startmol\" is not a valid selection. Possible choices are $availmods"
   return
  }
  if {[checkexists $startmol]} {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Starting model \"$startmol\" (or elements of it) is no longer in memory. Use the Refresh button to update list of available selections."
   return
  }
  if {[lsearch $availmods $endmol] == -1 } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Ending model \"$endmol\" is not a valid selection. Possible choices are $availmods"
   return
  }
  if {[checkexists $endmol]} {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Ending model \"$endmol\" (or elements of it) is no longer in memory. Use the Refresh button to update list of available selections."
   return
  }
 
  if {$endmol == $startmol } {
   tk_messageBox -icon error -type ok -title Message -parent $wan \
   -message "Same model can not be used for both the starting and ending point."
   return

  }
  menu tkcon on
  puts "\nWill generate an animation from molecule animation ID $startmol to $endmol."
  ::RTANIMATE::RT-animate2 $startmol $endmol $stepnum $MODN $VSEQ 
 }

 proc checkexists {i} {
  set RTHIST $::RADTOOL::radenv(RADtoolhistory)
  set entry [lindex $RTHIST $i]
  set LSUID [lindex $entry end-1]
  set SSUID [lindex $entry end]
  set mod1 [lindex $entry end-4]
  set mod2 [lindex $entry end-3]
  if { (([lsearch [molinfo list] $LSUID] != -1 &&  [lsearch [molinfo list] $mod1] != -1  ) || $LSUID == -100 ) && [lsearch [molinfo list] $SSUID] != -1 && [lsearch [molinfo list] $mod2] != -1 } {
   # still in memory
   return 0
  } else {
   # not valid
   return 1
  }
 }

 proc reset {} {
  set ::RADGUI::RCSBname {}
  set ::RADGUI::LSUSSU {all}
  set ::RADGUI::doSTAMP {1}
  set ::RADGUI::doDUMP {0}
  set ::RADGUI::pruneON {1}
  set ::RADGUI::pruneby {2}
  set ::RADGUI::animateON {none}
  set ::RADGUI::LSUfile {} 
  set ::RADGUI::SSUfile {}
  set ::RADGUI::LSUchain {}
  set ::RADGUI::SSUchain {}
  set ::RADGUI::saveDir {~/}
  set ::RADGUI::savePrefix {}
  set ::RADGUI::BUNDLEfile {}
  set ::RADGUI::BUNDLEfilestr ""
 }

 proc radmenu {} {
  RAD-GUI 
  return $::RADGUI::w
 }

}

