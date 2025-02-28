#pragma once

#include <cstdint>

void trs4slip_astar(
    int32_t * x_next_out,
    double const * c,
    int32_t const * x,
    int32_t const * bangs,
    int const delta,
    int const N,
    int const M,
    double * vert_costs_buffer,
    int32_t * vert_layer_buffer,
    int32_t * vert_value_buffer,
    int32_t * vert_prev_buffer,
    int32_t * vert_remcap_buffer,
    int boundcon,
    double lbound,
    double rbound 
);

void trs4slip_top(
    int32_t * x_next_out,
    double const * c,
    int32_t const * x,
    int32_t const * bangs,
    double const * switchingcost,
    int const Delta,
    int const N,
    int const M,
    int actbounds, 
    double leftbound,
    double rightbound,
    double switchcostleft,
    double switchcostright
);
