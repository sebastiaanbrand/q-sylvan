/*
 * Copyright 2023 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef QSYLVAN_MTBDD_GATES_H
#define QSYLVAN_MTBDD_GATES_H

#include <sylvan_mtbdd.h>

/**
 *
 * Definition of the gate matrices 2 x 2, with mpc floatingpoint type.
 *
 */

typedef enum quantum_gate {
    I,
    X,
    Y,
    Z,
    H,
    S,
    Sdag,
    T,
    Tdag,
    sqrtX,
    sqrtXdag,
    sqrtY,
    sqrtYdag
} quantum_gate_t;

MTBDD quantum_gate_to_mtbdd(quantum_gate_t g);

#endif
