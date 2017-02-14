/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016, 2017 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <iostream>
#include "data.h"
#include "vertical_segment.h"
#include "common.h"
#include "args.h"

namespace alta
{
    // Write DATA to OUT in ALTA's text format.
    void save_data_as_text(std::ostream& out, const alta::data &data);

    // Write DATA to OUT in a compact binary format.
    void save_data_as_binary(std::ostream& out, const alta::data& data);


    // Return the data read from INPUT in ALTA's native text format.  Read
    // extra parsing parameters from ARGS.  XXX: ARGS is currently used for
    // extra filtering, which should really happen elsewhere.
    data* load_data_from_text(std::istream& input,
                              const alta::arguments& header,
                              const alta::arguments& args = alta::arguments());

    // Return the data read from the binary-formatted stream IN.
    data* load_data_from_binary(std::istream& in, const alta::arguments& header);
}

