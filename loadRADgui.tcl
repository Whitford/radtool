# this is just a short wrapper to load RADtool using the GUI
set rtroot [ file dirname [ file normalize [ info script ] ] ]
source $rtroot/tcl/RADtool.tcl
source $rtroot/tcl/RAD-GUI.tcl
namespace import RADTOOL::RADtool
puts stderr "Launching RAD-GUI

If you would like to have RADtool persistently available through the VMD Menu, then add the following line to your .vmdrc or vmd.rc file (see README):

source \"$rtroot/installRADgui.tcl\"

Then close and reopen VMD. If you added this to the same rc file that is read by VMD, then the GUI will be available through the VMD Main menu item Extensions->Analysis->RADtool"

::RADGUI::RAD-GUI
