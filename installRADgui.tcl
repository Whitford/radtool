# this is just a short wrapper to load RADtool using the GUI
set rtroot [ file dirname [ file normalize [ info script ] ] ]
source $rtroot/tcl/RADtool.tcl
source $rtroot/tcl/RAD-GUI.tcl
namespace import RADTOOL::RADtool
lappend auto_path $rtroot/tcl
vmd_install_extension radgui ::RADGUI::radmenu  "Analysis/RADtool"
