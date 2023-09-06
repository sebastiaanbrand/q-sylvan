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

#include <qsylvan_mtbdd_gates.h>

MTBDD gate_to_mtbdd(quantum_gate_t g)
{
    if(g == I) {
        mpc_t G[2][2];
        G[0][0] = MTBDD_ONE;  G[0][1] = MTBDD_ZERO;
        G[1][0] = MTBDD_ZERO; G[1][1] = MTBDD_ONE;
        return array_matrix_to_mtbdd(G,2);
    }

    // ...

    return MTBDD_ZERO;
}

/*
    Pi = 2.0 * flt_acos(0.0);

    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_ONE;

    k = GATEID_X;
    gates[k][0] = AADD_ZERO; gates[k][1] = AADD_ONE;
    gates[k][2] = AADD_ONE;  gates[k][3] = AADD_ZERO;

    k = GATEID_Y;
    gates[k][0] = AADD_ZERO; gates[k][1] = weight_lookup(cmake(0.0, -1.0));
    gates[k][2] = weight_lookup(cmake(0.0, 1.0));  gates[k][3] = AADD_ZERO;

    k = GATEID_Z;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = weight_lookup(cmake(1.0/flt_sqrt(2.0),0));
    gates[k][3] = weight_lookup(cmake(-1.0/flt_sqrt(2.0),0));

    k = GATEID_S;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, 1.0));

    k = GATEID_Sdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, -1.0));

    k = GATEID_T;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5,-0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));

    k = GATEID_sqrtXdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5, 0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));

    k = GATEID_sqrtY;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(-0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));

    k = GATEID_sqrtYdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(-0.5,0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));

    qmdd_phase_gates_init(255);

    next_custom_id = 0;

*/

