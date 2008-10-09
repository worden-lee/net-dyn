/* -*- c++ -*-
 */
#ifndef LPARAMETERS_H
#define LPARAMETERS_H

#include "../network/Parameters.h"

#define BLOCKSIZE 128
#define NBLOCKS 8

#undef DISPLAY

class LParameters : public Parameters
{
public:
  DECLARE_PARAM(string, desKey);
};

#endif//LPARAMETERS_H
