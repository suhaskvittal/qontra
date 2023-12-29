/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef LINPROG_MANAGER_h
#define LINPROG_MANAGER_h

#include <cplexx.h>

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
    bool solve(double* objective_p, int* solution_status_p);

    void solve_pool(void);

    lp_var_t    get_var(T);
    double      get_value(T);

    lp_var_t    add_slack_var(double lower_bound, double upper_bound, lp_var_t::bounds, lp_var_t::domain);
    lp_var_t    add_var(T label, double lower_bound, double upper_bound, lp_var_t::bounds, lp_var_t::domain);

    size_t  add_constraint(lp_constr_t);

    size_t  get_soln_pool_size(void);
    double  fetch_soln_from_pool(size_t);

    std::vector<lp_var_t>    variables;
    std::vector<lp_constr_t> constraints;
private:
    enum class problem_type { lp, mip };

    std::map<T, lp_var_t> label_to_lp_var;

    size_t columns;
    size_t rows;

    CPXENVptr env;
    CPXLPptr prog;
    bool env_is_initialized_by_object;

    problem_type prog_type;

    double* prog_soln;
};

#include "manager.inl"

#endif  // LINPROG_MANAGER_h
