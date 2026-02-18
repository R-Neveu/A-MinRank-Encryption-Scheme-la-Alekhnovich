/**
 * \file lrpc.h
 * \brief Fonctions to decode the error support using LRPC codes
 *
 */

#ifndef RBC_61_LRPC_H
#define RBC_61_LRPC_H

#include "rbc_61_vspace.h"
#include "rbc_61_qre.h"

unsigned int rbc_61_lrpc_RSR(rbc_61_vspace E, const uint32_t E_expected_dim, const rbc_61_vspace F, const uint32_t F_dim, const rbc_61_vec ec, const uint32_t ec_size);

#endif

