#!/usr/bin/env python
# $Id: LQanalyzer_create_cycle.py jalmond $
#***************************************************************************
#* @Project: LQanalyzer - ROOT-based analysis framework for CMS
#* @Package: Core
#*
#* @author J. Almond  <jalmond@cern.ch> - SNU
#*
#***************************************************************************
#
# This script can be used to quickly create a new analysis cycle in
# the user's package. If it's called in an arbitrary directory, then
# all the created files are put in the current directory. However
# when invoking it in an LQanalzer library directory (like LQanalzer/user)
# it will put the files in the correct places and adds the entry
# to the already existing LinkDef file.
#

# Import base module(s):
import sys
import os.path
import optparse

## Main function of this script
#
# The executable code of the script is put into this function, much like
# it is done in the case of C++ programs. The main function is then called
# from outside using the usual method. (Which I've seen in various
# places...)
def main():
      # Print some welcome message before doing anything else:
      print ">>"
      print ">> %s : Analysis cycle torso creator" % \
            os.path.basename( sys.argv[ 0 ] )
      print ">>"
      print ""

      # Parse the command line parameters:
      parser = optparse.OptionParser( usage="%prog [options] <input files>" )
      parser.add_option( "-n", "--name", dest="name", action="store",
                         type="string", default="AnalysisCycle",
                         help="Name of the analysis cycle to create" )
      parser.add_option( "-l", "--linkdef", dest="linkdef", action="store",
                         type="string", default="",
                         help="Name of the LinkDef.h file in the package" )

      ( options, args ) = parser.parse_args()

      # This is where the main function are:
      import CycleCreators
      # Execute the cycle creation:
      cc = CycleCreators.CycleCreator()
      cc.CreateCycle( options.name, options.linkdef )

# Call the main function:
if __name__ == "__main__":
    main()
