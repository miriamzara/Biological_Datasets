#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)
# 2018.01.16 - Martin Voegele - Update to use Gromacs version 5 and later 

from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess

### SOME TOOLS

# parse a .gro file
# return a list of coordinates
def read_gro(file, atoms):
  line_counter = 0
  number_of_particles = 0
  first, second = [], []
  for line in open(file):
    if line_counter == 1:
      number_of_particles = int(line)
    elif line_counter > 1 and line_counter < number_of_particles + 2:
      if line[10:15].strip() == atoms[0]:
        first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
      elif line[10:15].strip() == atoms[1]:
        second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    line_counter += 1
  return [first, second]

### REAL STUFF

if len(argv) != 10:
  # coments/usage
  print("""
  Compute (second rank) order parameter, defined as:

    P2 = 0.5*(3*<cos²(theta)> - 1)

  where "theta" is the angle between the bond and the bilayer normal.
  P2 = 1      perfect alignement with the bilayer normal
  P2 = -0.5   anti-alignement
  P2 = 0      random orientation

  All lipids defined in the "martini_v2.0_lipids.itp" file can be analyzed
  with this script.
  Usage: %s <traj file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <#lipids> <lipid type>

    > %s traj.xtc 0 10000 5 0 1 0 64 DSPC

  will for example read a 10ns trajectory of 64 DSPC lipids, calculating the order parameter for 
  every 5th frame and averaging the results. P2 will be calculated relative to the y-axis.

  WARNING script will output all frames in one go, into files called frame_dump_XXX.gro and 
  then remove them so don't have any other files with this name in the current directory.
  """) % (argv[0], argv[0])
  exit(0)

else:

  # snapshots
  trajfile = argv[1]
  initial_time = int(argv[2])
  final_time = int(argv[3])
  traj_skip = int(argv[4])
  # (normalized) orientation of bilayer normal
  orientation_of_bilayer_normal = [float(argv[5]), float(argv[6]), float(argv[7])]
  norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
  for i in range(3):
    orientation_of_bilayer_normal[i] /= norm
  stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
    orientation_of_bilayer_normal[0], \
    orientation_of_bilayer_normal[1], \
    orientation_of_bilayer_normal[2]  \
  ))
  # number of lipids
  number_of_lipids = int(argv[8])
  # lipid type
  lipid_type = argv[9]

  # output legend
  phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
  # PCs
  if   lipid_type == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DHPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
  elif lipid_type == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DSPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  # PEs
  elif lipid_type == "DHPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
  elif lipid_type == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DSPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  # PPCS
  elif lipid_type == "PPCS": bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 AM1-C1A GL2-D1B C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B\n"
  # output legend
  output_legend = "  Frame" + bond_names 

  # write the stuff
  stdout.write("\n " + output_legend)
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")
  output = open('order.dat', 'w')
  output.write(output_legend)
  #output.write(("-"*(len(output_legend) - 1)) + "\n")

  # Output all frame using trjconv 
  stdout.write("Output all coordinate files \n")
  command = "echo %s | gmx trjconv -f %s -b %i -e %i -sep -skip %i -pbc whole -o frame_dump_.gro > /dev/null" % (lipid_type, trajfile, initial_time, final_time, traj_skip)
  print(command)
  subprocess.call(command, shell=True)

  # For each dumped frame
  stdout.write("Starting P2 calculation")
  order_parameters = []
  file_count = 0
  bonds = []
  while True:
    filename = "frame_dump_" + str(file_count) + ".gro"
    if not path.isfile(filename) or path.getsize(filename) == 0:
        break
    
    stdout.write("Taking care of snapshot %s \n" % filename)

    # compute order parameter for each bond, for each snapshot
    current_order_parameters = []
    # bonds respectively involved in the head,
    #                             in the junction head-tail,
    #                             in each tail
    bonds = []

    for bond_name in bond_names.split():
      bonds.append(bond_name.split("-"))

    for bond in bonds:

      # parse .gro file, grep bead coordinates
      first, second = read_gro(filename, bond)

      # compute order parameter for each lipid
      order_parameter = 0.0
      for i in range(number_of_lipids):
        # vector between the two previous beads (orientation doesn't matter)
        vector = [0.0, 0.0, 0.0]
        for j in range(3):
          vector[j] = first[i][j] - second[i][j]
        norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
        # compute projection on the bilayer normal
        projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
        # order parameter
        order_parameter += projection**2/norm2

      # compute final averaged order parameter
      # store everything in lists
      current_order_parameters.append(0.5*(3.0*(order_parameter/number_of_lipids) - 1.0))
    order_parameters.append(current_order_parameters)

    # write results
    results = "%7i" % file_count
    for order_parameter in current_order_parameters:
      results += "%8.3f" % order_parameter
    stdout.write(" " + results + "\n")
    output.write(results + "\n")

    remove(filename)
    file_count += 1
  # End while loop

  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  stdout.write("Snapshots analysis done.%s\n" % (" "*56))
  stdout.write("Computing averages...\n")

  # average order parameter
  averaged_order_parameters = []
  for i in range(len(bonds)):
    sum = 0.0
    for j in range(len(order_parameters)):
      sum += order_parameters[j][i]
    averaged_order_parameters.append(sum/len(order_parameters))
 
  # write results
  stdout.write("\n         " + bond_names)
  stdout.write(("-"*(len(output_legend) - 1)) + "\n")
  output.write(("-"*(len(output_legend) - 1)) + "\n")
  results = "average"
  for order_parameter in averaged_order_parameters:
    results += "%8.3f" % order_parameter
  stdout.write(" " + results + "\n")
  output.write(results + "\n")
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")

  # Write abs average order parameters <Sn> (for carbon chains only)
  # WARNING this works with currenct lipids (all have defined x5 none carbon bonds) but for manually added lipids this might not be true
  ave_chain_s = 0
  for i in averaged_order_parameters[3:]: 
     ave_chain_s += abs(i)
  average_txt = "Abs average order parameters for carbon chains <Sn> = %8.3f \n\n" % (ave_chain_s / (len(averaged_order_parameters)-3))
  stdout.write(average_txt)
  output.write(average_txt)
  stdout.write("Results written in \"order.dat\".\n")

  output.close()
