/**
 * \file lrpc.h
 * \brief Fonctions to decode the error support using LRPC codes
 *
 */

#ifndef RBC_37_LRPC_H
#define RBC_37_LRPC_H

#include "rbc_37_vspace.h"
#include "rbc_37_qre.h"

unsigned int rbc_37_lrpc_RSR(rbc_37_vspace E, const uint32_t E_expected_dim, const rbc_37_vspace F, const uint32_t F_dim, const rbc_37_vec ec, const uint32_t ec_size);

#endif

