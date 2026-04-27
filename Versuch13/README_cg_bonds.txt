Use cg_bonds.tcl to visualize MARTINI bonds in VMD.

Find the path to program 'gmxdump', execute 

$ which gmxdump

in a shell.

Start vmd and type

source cg_bonds.tcl
cg_bonds -gmx PATH -tpr TOPO

where PATH is the path to gmxdump and TOPO is the name of your topolgy file (.tpr).  
