
#ifndef VTILS_LINPROG_MANAGER_h
#define VTILS_LINPROG_MANAGER_h

#include "vtils/linprog/base.h"

#include <cplexx.h>

namespace vtils {

void handle_status(int);

CPXENVptr   cpxinit(void);
void        cpxexit(CPXENVptr*);  
void        cpxmakeprog(CPXENVptr, CPXLPptr, const char* name);
void        cpxfreeprog(CPXENVptr, CPXLPptr);

template <class T>
class CPXLPManager {
public:
    CPXLPManager(CPXENVptr env=NULL);
    ~CPXLPManager(void);

    void build(lp_expr_t objective, bool is_maximization);
    // Returns nonzero if solve fails.
    int solve(double* objective_p, int* solution_status_p);

    void solve_pool(void);

    lp_var_t    get_var(T);
    double      get_value(lp_var_t);
    double      get_value(T);

    lp_var_t    add_slack_var(double lower_bound, double upper_bound, lp_var_t::bounds, lp_var_t::domain);
    lp_var_t    add_var(T label, double lower_bound, double upper_bound, lp_var_t::bounds, lp_var_t::domain);

    void    add_constraint(lp_constr_t);

    size_t  get_soln_pool_size(void);
    double  fetch_soln_from_pool(size_t);

    std::vector<lp_var_t>    variables;
    std::vector<lp_constr_t> l_constraints;
    std::vector<lp_constr_t> q_constraints;
private:
    char get_sense(lp_constr_t);

    enum class problem_type { lp, mip, qp, qcp};

    std::map<T, lp_var_t> label_to_lp_var;

    CPXENVptr env;
    CPXLPptr prog;

    problem_type prog_type;

    double* prog_soln;

    bool env_is_initialized_by_object;
};

}   // vtils

#include "inl/manager.inl"

#endif  // LINPROG_MANAGER_h
