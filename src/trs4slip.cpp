#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <queue>

using namespace std;

#include "trs4slip.h"

struct vertex_rel_t
{
    double cost;
    int used;
    int prev;
};

typedef vector<vector<vertex_rel_t>> cost_dict_t;

struct vertex_t
{
    double cost;
    int idx;
};

double const BIN_SEARCH_EPS = 1e-9; 
int const NUM_RUNS_MAX = 30;        // depends on BIN_SEARCH_EPS
                                    // and be reduced if BIN_SEARCH_EPS is smaller

bool binsearch(
    int32_t * x_next_out,
    cost_dict_t & cost_dict_out,
    vector<double> & penaltylist_out,
    int & penalty_num,
    int & num_runs,
    double & ub,
    int delta,
    int N,
    int M,
    double const *c,
    int32_t const *x,
    int32_t const *bangs,
    double offset,
    int boundcon,
    double lbound,
    double rbound
);

inline double addcost(
    int N,
    int d,
    int pred,
    int32_t const *x,
    double const *c,
    int layer,
    int boundcon,
    double lbound,
    double rbound)
{
    if (layer > 1)
    {
        if (layer == N and boundcon)
        {
            return c[layer - 1] * d + abs(x[layer - 1] + d - pred) + abs(x[layer - 1] + d - rbound);
        }
        return c[layer - 1] * d + abs(x[layer - 1] + d - pred);
    }
    if (boundcon)
    {
        return c[layer - 1] * d + abs(x[layer - 1] + d - lbound);
    }
    return c[layer - 1] * d;
}

void trs4slip_astar(
    int32_t * x_next_out,
    double const *c,
    int32_t const *x,
    int32_t const *bangs,
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
)
{
    assert(delta >= 0.);

    double offset = abs(c[N-1]) * max(x[N-1] - bangs[0], bangs[M-1] - x[N-1]);
    for (int i = N - 2; i >= 0; i--)
    {
        double const max_element = abs(c[i]) * max(x[i] - bangs[0], bangs[M-1] - x[i]);
        if (offset < max_element)
            offset = max_element;
    }

    vector<double> tv(N);
    tv[N - 1] = 0.;
    if (boundcon)
    {
        tv[N - 1] = abs(rbound - x[N - 1]);
    }
    for (int i = 1; i < N; i++)
        tv[N - i - 1] = abs(x[N - i] - x[N - i - 1]) + tv[N - i];

    vector<double> penaltylist;
    penaltylist.reserve(NUM_RUNS_MAX);

    cost_dict_t cost_dict(N * M + 1);
    for (int i = 0; i < N * M + 1; i++)
        cost_dict.reserve(NUM_RUNS_MAX);

    int penalty_num = 0;
    int num_runs = NUM_RUNS_MAX;
    double ub = numeric_limits<double>::max();

    if (binsearch(&x_next_out[0], cost_dict, penaltylist, penalty_num, num_runs, ub,
                  delta, N, M, &c[0], &x[0], &bangs[0], offset, boundcon, lbound, rbound))
        return;

    fill(&vert_costs_buffer[0], &vert_costs_buffer[N * M * (delta + 1) + 2], numeric_limits<double>::max());
    std::vector<bool> visited(N * M * (delta + 1) + 2, false);

    double penalty;
    int num_binsearch = 0;
    int num_astar = 0;
    vert_costs_buffer[0] = 0.0;
    vert_layer_buffer[0] = 0;
    vert_prev_buffer[0] = 0;
    vert_value_buffer[0] = 0;
    vert_remcap_buffer[0] = delta;

    auto compare = [](vertex_t const &a, vertex_t const &b)
    { return a.cost > b.cost; };
    priority_queue<vertex_t, vector<vertex_t>, decltype(compare)> prio(compare);
    prio.push({.cost = 0.0, .idx = 0});

    double max_lowerbound;
    double trivial_achievable_cost;

    while (!prio.empty())
    {
        int current = prio.top().idx;
        prio.pop();
        while (visited[current] and !prio.empty())
        {
            current = prio.top().idx;
            prio.pop();
        }
        visited[current] = true;
        int layer = vert_layer_buffer[current];
        int const curr_remcap = vert_remcap_buffer[current];
        int const curr_value = vert_value_buffer[current];
        if (layer == N)
        {
            int idx = current;
            for (int i = N; i > 0; i--)
            {
                x_next_out[i - 1] = bangs[vert_value_buffer[idx]];
                idx = vert_prev_buffer[idx];
            }
            return;
        }
        layer++;
        for (int k = 0; k < M; k++)
        {
            int const d = bangs[k] - x[layer - 1];
            bool improve = false;
            int const remcap = curr_remcap - abs(d);
            if (remcap >= 0)
            {
                double const cost = vert_costs_buffer[current] + offset + addcost(N, d, bangs[curr_value], x, c, layer, boundcon, lbound, rbound);
                int const vert_idx_in_next_layer = ((layer - 1) * (delta + 1) + remcap) * M + k + 1;
                if (vert_costs_buffer[vert_idx_in_next_layer] > cost)
                {
                    vert_costs_buffer[vert_idx_in_next_layer] = cost;
                    max_lowerbound = 0.0;
                    if (layer<N and layer > 1 and remcap > 0)
                    {
                        int const vert_idx_in_binsearch = (layer - 1) * M + k + 1;
                        for (int t = num_runs - 1; t >= 0; t--)
                        {
                            penalty = penaltylist[t];
                            vertex_rel_t const &cost_ref = cost_dict[vert_idx_in_binsearch][t];
                            max_lowerbound = max(max_lowerbound, -penalty * remcap + cost_ref.cost);
                            if (cost_ref.used <= remcap)
                            {
                                if (ub > cost_ref.cost - penalty * cost_ref.used + cost)
                                {
                                    ub = cost_ref.cost - penalty * cost_ref.used + cost;
                                    penalty_num = t;
                                    num_binsearch = vert_idx_in_binsearch;
                                    num_astar = vert_idx_in_next_layer;
                                    vert_layer_buffer[vert_idx_in_next_layer] = layer;
                                    vert_value_buffer[vert_idx_in_next_layer] = k;
                                    vert_remcap_buffer[vert_idx_in_next_layer] = remcap;
                                    vert_prev_buffer[vert_idx_in_next_layer] = current;
                                }
                            }
                            if (max_lowerbound + cost > ub + max(1e-12, 1e-5 * abs(ub)) )
                            {   

                                improve = true;
                                break;
                            }
                        }
                    }
                    if (remcap == 0 and layer < N)
                    {
                        trivial_achievable_cost =
                            cost + tv[layer] + abs(x[layer] - x[layer - 1] - d) + (N - layer) * offset;
                        if (trivial_achievable_cost > ub + max(1e-12, 1e-5 * abs(ub)))
                            improve = true;
                        if (ub > trivial_achievable_cost)
                        {
                            ub = trivial_achievable_cost;
                            penalty_num = -1;
                            num_binsearch = (layer - 1) * M + k + 1;
                            num_astar = vert_idx_in_next_layer;
                            vert_layer_buffer[vert_idx_in_next_layer] = layer;
                            vert_value_buffer[vert_idx_in_next_layer] = k;
                            vert_remcap_buffer[vert_idx_in_next_layer] = remcap;
                            vert_prev_buffer[vert_idx_in_next_layer] = current;
                        }
                    }
                    if (!improve)
                    {
                        vert_layer_buffer[vert_idx_in_next_layer] = layer;
                        vert_value_buffer[vert_idx_in_next_layer] = k;
                        vert_remcap_buffer[vert_idx_in_next_layer] = remcap;
                        vert_prev_buffer[vert_idx_in_next_layer] = current;
                        if (layer == N)
                            prio.push({.cost = cost, .idx = vert_idx_in_next_layer});
                        else
                        {
                            if (remcap == 0)
                                prio.push({.cost = cost + tv[layer] + abs(x[layer] - x[layer - 1] - d) + (N - layer) * offset, .idx = vert_idx_in_next_layer});
                            else
                                prio.push({.cost = cost + max_lowerbound, .idx = vert_idx_in_next_layer});
                        }
                    }
                }
            }
            else if (d > 0)
                break;
        }
    }

    std::cout << "Warning: Rounding Errors, astar failed, returned solution might be suboptimal" << std::endl;

    int iter_num = num_astar;
    int const vert_layer_start = vert_layer_buffer[iter_num];
    int vl = vert_layer_start;
    for (; iter_num > 0; iter_num = vert_prev_buffer[iter_num])
        x_next_out[--vl] = bangs[vert_value_buffer[iter_num]];

    if (penalty_num < 0)
        copy(&x[vert_layer_start], &x[N], &x_next_out[vert_layer_start]);
    else
    {
        vl = vert_layer_start;
        if (num_binsearch > 0 and num_astar <= 0)
        {
            x_next_out[0] = bangs[(iter_num - 1) % M];
            vl = 1;
        }
        iter_num = num_binsearch;
        for (; vl < N; vl++)
        {
            iter_num = cost_dict[iter_num][penalty_num].prev;
            x_next_out[vl] = bangs[(iter_num - 1) % M];
        }
    }
}

bool binsearch(
    int32_t * x_next_out,
    cost_dict_t & cost_dict_out,
    vector<double> & penaltylist_out,
    int & penalty_num, 
    int & num_runs,
    double & ub,
    int delta,
    int N,
    int M,
    double const *c,
    int32_t const *x,
    int32_t const *bangs,
    double offset,
    int boundcon,
    double lbound,
    double rbound)
{
    vector<double> cmax(N);
    cmax[N - 1] = abs(c[N - 1]);
    for (int i = N - 2; i >= 0; i--)
    {
        double const max_element = abs(c[i]);
        if (cmax[i + 1] > max_element)
            cmax[i] = cmax[i + 1];
        else
            cmax[i] = max_element;
    }

    vector<double> relaxed_costs(N * M + 1);
    vector<int> relaxed_prev(N * M + 1);
    vector<int> relaxed_used(N * M + 1);
    vector<int> relaxed_value(N * M + 1);

    assert(cmax[0] >= 0.);

    double upper = cmax[0] + 2. + 1e-5;
    double lower = 0.0;
    bool optimal = false;
    for (int t = 0; t < num_runs; t++)
    {
        if (upper - lower < BIN_SEARCH_EPS)
        {
            num_runs = t;
            return false;
        }

        double const penalty = (upper + lower) / 2.0;
        for (int k = 0; k < M; k++)
        {
            relaxed_costs[N * M - k] = 0;
            if (boundcon)
            {
                relaxed_costs[N * M - k] = abs(rbound - bangs[M - 1 - k]);
            }
            relaxed_prev[N * M - k] = -1;
            relaxed_used[N * M - k] = 0;
            relaxed_value[N * M - k] = M - 1 - k;
        }

        for (int k = 0; k < N - 1; k++)
        {
            for (int l = 0; l < M; l++)
            {
                int const cost_dict_idx = (N - k) * M - l;
                cost_dict_out[cost_dict_idx].push_back(
                    {.cost = relaxed_costs[cost_dict_idx],
                     .used = relaxed_used[cost_dict_idx],
                     .prev = relaxed_prev[cost_dict_idx]});
                if (l == 0)
                {
                    int const d = bangs[M - 1 - l] - x[N - 1 - k];
                    for (int o = 0; o < M; o++)
                    {
                        int const im_cost_idx = (N - 2 - k) * M + o + 1;
                        relaxed_costs[im_cost_idx] =
                            relaxed_costs[cost_dict_idx] + abs(bangs[M - 1 - l] - bangs[o]) + c[N - 1 - k] * d + offset + abs(d) * penalty;
                        relaxed_prev[im_cost_idx] = cost_dict_idx;
                        relaxed_used[im_cost_idx] = relaxed_used[cost_dict_idx] + abs(d);
                        relaxed_value[im_cost_idx] = o;
                    }
                }
                else
                {
                    int const d = bangs[M - 1 - l] - x[N - 1 - k];
                    for (int o = 0; o < M; o++)
                    {
                        int const im_cost_idx = (N - 2 - k) * M + o + 1;
                        double const cost =
                            relaxed_costs[cost_dict_idx] + abs(bangs[M - 1 - l] - bangs[o]) + c[N - 1 - k] * d + offset + abs(d) * penalty;
                        if (cost < relaxed_costs[im_cost_idx])
                        {
                            relaxed_costs[im_cost_idx] = cost;
                            relaxed_prev[im_cost_idx] = (N - k) * M - l;
                            relaxed_used[im_cost_idx] = relaxed_used[cost_dict_idx] + abs(d);
                            relaxed_value[im_cost_idx] = o;
                        }
                    }
                }
            }
        }
        for (int k = 0; k < M; k++)
        {
            relaxed_costs[k + 1] =
                relaxed_costs[k + 1] + c[0] * (bangs[k] - x[0]) + penalty * abs(bangs[k] - x[0]) + offset;
            relaxed_used[k + 1] = relaxed_used[k + 1] + abs((bangs[k] - x[0]));
            cost_dict_out[k + 1].push_back(
                {.cost = relaxed_costs[k + 1],
                 .used = relaxed_used[k + 1],
                 .prev = relaxed_prev[k + 1]});
        }
        int max_cost_idx = 1;
        double max_cost = relaxed_costs[max_cost_idx] ;
        if (!boundcon) {
            
            for (int k = 1; k < M; k++)
            {
                if (max_cost > relaxed_costs[k + 1] )
                {
                    max_cost = relaxed_costs[k + 1] ;
                    max_cost_idx = k + 1;
                }
            }
        }
        if (boundcon) {
            max_cost_idx = 1;
            max_cost = relaxed_costs[max_cost_idx] + abs(bangs[0] -lbound);
            for (int k = 1; k < M; k++)
            {
                if (max_cost > relaxed_costs[k + 1] + abs(bangs[k] - lbound ))
                {
                    max_cost = relaxed_costs[k + 1] + abs(bangs[k] - lbound );
                    max_cost_idx = k + 1;
                }
            }
        }
        cost_dict_out[0].push_back(
            {.cost = max_cost,
             .used = relaxed_used[max_cost_idx],
             .prev = max_cost_idx});

        penaltylist_out.push_back(penalty);

        if (relaxed_used[max_cost_idx] > delta)
            lower = penalty;
        else
        {
            if (relaxed_used[max_cost_idx] == delta)
            {
                optimal = true;
                for (int ii = 0; ii < N; ii++)
                {
                    x_next_out[ii] = bangs[relaxed_value[max_cost_idx]];
                    max_cost_idx = relaxed_prev[max_cost_idx];
                }
                return optimal;
            }
            if (ub > max_cost - relaxed_used[max_cost_idx] * penalty)
            {
                ub = max_cost- relaxed_used[max_cost_idx] * penalty;
                penalty_num = t;
            }
            upper = penalty;
        }
    }
    return optimal;
}

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
    double switchcostright) {

    int d;
    int usedcap;
    int index;
    int predindex;
    int dabs;
    double bestcost;
    int Deltaplus1 = Delta+1;
    int nmalm = N*M;
    double finalcost =  std::numeric_limits<double>::max();
    int bestvert=0;
    double nodecost;
    double newcost;
    double precost;
    vector<double> costvec1(nmalm*(Deltaplus1), std::numeric_limits<double>::max());
    vector<int> prevvec1(nmalm*(Deltaplus1), -1);
    vector<int> valuevec1(nmalm*(Deltaplus1), 0);   
    
    for (int o = 0; o < M; o++) {
        d = bangs[o] - x[0];
        usedcap = abs(d);
        if (usedcap <= Delta) {
            index = o * Deltaplus1 + usedcap;
            if (actbounds ==1) {
                costvec1[index] = c[0] * d+ switchcostleft*abs(x[0] + d - leftbound);
            } else {
                costvec1[index] = c[0] * d;
            }
            
            valuevec1[index] = o;
        }
    }
    for (int k = 1; k < N; k++){
        int prevk = k-1;
        for (int o = 0; o < M; o++) {
            for (int p = 0; p < M; p++) {
                d = bangs[p] - x[k];
                precost = c[k]*d;
                precost = precost + switchingcost[prevk] * abs(bangs[p]-bangs[o]);
                if (k==N-1 && actbounds== 1) {
                    precost = precost + switchcostright*abs(bangs[p]- rightbound);
                }
                dabs = abs(bangs[p]-x[k]);
                index = (k*M+p)*Deltaplus1+dabs;
                
                predindex = (prevk*M+o)*Deltaplus1;
                bestcost = std::numeric_limits<double>::max();
                for (int r = Delta-dabs; r >= 0; r--) {
                    nodecost = costvec1[predindex];
                    if (nodecost<bestcost) {
                        //bestcost = nodecost;
                        newcost = nodecost + precost;
                        if ( newcost<costvec1[index]) {
                            costvec1[index] = newcost;
                            prevvec1[index] = predindex;
                            valuevec1[index] = p;
                            if (k==N-1) {
                                if (newcost<finalcost) {     
                                    finalcost = newcost;
                                    bestvert = index;
                                }
                            }
                        }
                    }
                    predindex++;
                    index++;
                }
            }
            
        }
    }

    int p;
    int ctr = 0;
    for (int j=N-1; j>=0; j--) {
        p = valuevec1[bestvert];
        x_next_out[ctr++] = bangs[p];
        bestvert = prevvec1[bestvert];
    }
	std::reverse(&x_next_out[0], &x_next_out[N]);
}
