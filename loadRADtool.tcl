# this is just a short wrapper to load RADtool without the GUI
set rtroot [ file dirname [ file normalize [ info script ] ] ]
source $rtroot/tcl/RADtool.tcl
namespace import RADTOOL::RADtool
