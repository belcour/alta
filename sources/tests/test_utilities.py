# ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

#   Copyright (C) 2017, CNRS

#   This file is part of ALTA.

#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0.  If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#   

#   author:  Romain Pacanowski:  romain.pacanowski@institutoptique.fr
#  



##############################################################################################
# Construct a dictionary  that contains the following keys (string type)
#  L1, LInf,
# TODO: CHECK FOR CONVERSION ERROR ?
# ##############################################################################################
def read_brdf_fit_file( brdf_fit_file ):
  file_handle = open(brdf_fit_file,  'r')
  row = file_handle.readlines();

  fit_info = {}

  for line in row:
    #print 'CURRENT:' + line + "____"
    #print 'Current print_info' + str(fit_info)

    if line.startswith("#CMD"):
      #TODO: get the L2 and LINF norm
      tokens = line.split(' ')

      if '--L2' in tokens:
        l2_index = tokens.index('--L2')
        fit_info['--L2'] =  float(tokens[l2_index+1])
      #end if

      if '--Linf' in tokens:
        linf_index = tokens.index('--Linf')
        fit_info['--Linf'] = float(tokens[linf_index+1])
      #end if

    elif line.startswith("#FUNC"):
      func_token, func_name = line.split(' ')
      fit_info[func_name.strip('\n')] = True
    elif (line != "\n" and not(line.startswith("#")) ):


      token, token_value = line.split()

      if fit_info.has_key(token):
        fit_info[token] += [float(token_value.strip('\n'))]
      else:
        fit_info[token] = [float(token_value.strip('\n'))]
      #end_if_else
    #end_if_else
  #end_for

  return fit_info
#end_def