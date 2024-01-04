/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#include <iostream>

#include <stdlib.h>

inline void handle_status(int s) {}

inline CPXENVptr cpxinit() {
    int status;
    CPXENVptr env = CPXXopenCPLEX(&status);
    handle_status(status);
    return env;
}

inline void cpxexit(CPXENVptr* env_p) {
    CPXXcloseCPLEX(env_p);
}

inline void cpxmakeprog(CPXENVptr env, CPXLPptr* prog_p, const char* name) {
    int status;
    *prog_p = CPXXcreateprob(env, &status, name);
    handle_status(status);
}

inline void cpxfreeprog(CPXENVptr env, CPXLPptr* prog_p) {
    if (*prog_p != NULL) {
        int status = CPXXfreeprob(env, prog_p);
        *prog_p = NULL;
        handle_status(status);
    }
}

template <class T>
CPXLPManager<T>::CPXLPManager(CPXENVptr _env)
    :variables(),
    l_constraints(),
    q_constraints(),
    label_to_lp_var(),
    env(_env),
    prog(NULL),
    prog_soln(NULL),
    env_is_initialized_by_object(_env == NULL)
{
    if (env_is_initialized_by_object) {
        env = cpxinit();

        CPXXsetdblparam(env, CPXPARAM_MIP_Pool_AbsGap, 0.0);
        CPXXsetintparam(env, CPXPARAM_MIP_Pool_Intensity, 4);
        CPXXsetintparam(env, CPXPARAM_MIP_Limits_Populate, 32);
    }
}

template <class T>
CPXLPManager<T>::~CPXLPManager() {
    cpxfreeprog(env, &prog);
    if (env_is_initialized_by_object) cpxexit(&env);
}

template <class T> void
CPXLPManager<T>::build(lp_expr_t objective, bool is_maximization) {
    int status;

    cpxfreeprog(env, &prog);
    cpxmakeprog(env, &prog, "lp");

    status = CPXXchgobjsen(env, prog, is_maximization ? CPX_MAX : CPX_MIN);
    handle_status(status);

    // Declare variables column-wise.
    const size_t columns = variables.size();
    double* obj = reinterpret_cast<double*>(calloc(columns, sizeof(double)));
    double* lb = reinterpret_cast<double*>(calloc(columns, sizeof(double)));
    double* ub = reinterpret_cast<double*>(calloc(columns, sizeof(double)));
    char* vtypes = reinterpret_cast<char*>(calloc(columns, sizeof(char)));

    // We should also track the problem type when constructing the problem.
    prog_type = objective.is_quadratic() ? problem_type::qp : problem_type::lp;
    if (q_constraints.size()) prog_type = problem_type::qcp;
    for (size_t i = 0; i < columns; i++) {
        lp_var_t v = variables[i];
        if (objective.l_coefs.count(v)) {
            obj[i] = objective.l_coefs[v];
        }

        bool lwr_defined = v.bounds_type == lp_var_t::bounds::lower
                            || v.bounds_type == lp_var_t::bounds::both
                            || v.bounds_type == lp_var_t::bounds::fixed;
        bool upp_defined = v.bounds_type == lp_var_t::bounds::upper
                            || v.bounds_type == lp_var_t::bounds::both
                            || v.bounds_type == lp_var_t::bounds::fixed;
        
        lb[i] = lwr_defined ? v.lwr : -CPX_INFBOUND;
        ub[i] = upp_defined ? v.upp : CPX_INFBOUND;
        if (v.var_type == lp_var_t::domain::continuous) {
            vtypes[i] = 'C';
        } else if (v.var_type == lp_var_t::domain::integer) {
            vtypes[i] = 'I';
            prog_type = problem_type::mip;
        } else {
            vtypes[i] = 'B';
            prog_type = problem_type::mip;
        }
    }
    status = CPXXnewcols(env, prog, columns, obj, lb, ub, vtypes, NULL);
    handle_status(status);
    // Free all allocated memory.
    free(obj);
    free(lb);
    free(ub);
    free(vtypes);
    // Set any quadratic coefficients if necessary.
    if (objective.is_quadratic()) {
        for (auto& pair : objective.q_coefs) {
            lp_quad_var_t qvar = pair.first;
            double x = pair.second;
            CPXXchgqpcoef(env, prog, qvar.v1.column, qvar.v2.column, x);
        }
    }

    // Now, declare the constraints.
    //
    // First, we need to count the total number of nonzero coefficients in 
    // the LINEAR constraints.
    int num_nonzeros = 0;
    for (auto& con : l_constraints) {
        num_nonzeros += con.lhs.l_coefs.size();
    }
    const size_t l_rows = l_constraints.size();
    // Now, we need to populate the data structures to declare the
    // constraints.
    double* rhs = reinterpret_cast<double*>(malloc(l_rows * sizeof(double)));
    // GE ('G'), LE ('L'), or EQ ('E').
    char* sense = reinterpret_cast<char*>(malloc(l_rows * sizeof(char)));
    // Where the coefficients for the i-th row.
    CPXNNZ* rmatbeg = reinterpret_cast<CPXNNZ*>(malloc(l_rows * sizeof(CPXNNZ)));
    // The column of the nonzero coefficient.
    int* rmatind = reinterpret_cast<int*>(malloc(num_nonzeros * sizeof(int)));
    // The value of the corresponding nonzero coefficient.
    double* rmatval = reinterpret_cast<double*>(malloc(num_nonzeros * sizeof(double)));

    size_t offset = 0;
    for (size_t i = 0; i < l_rows; i++) {
        lp_constr_t con = l_constraints[i];
        rhs[i] = con.rhs;
        sense[i] = get_sense(con);
        // We will not have anything else.
        rmatbeg[i] = offset;
        for (auto kv : con.lhs.l_coefs) {
            lp_var_t v = kv.first;
            rmatind[offset] = v.column;
            rmatval[offset] = kv.second;
            offset++;
        }
    }
    status = CPXXaddrows(env, prog, 0, l_rows, num_nonzeros, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
    handle_status(status);
    // Free all variables.
    free(rhs);
    free(sense);
    free(rmatbeg);
    free(rmatind);
    free(rmatval);
    // Now, if there are any quadratic constraints, we handle that now:
    //
    // Unlike linear constraints, these are handled one-by-one.
    for (lp_constr_t& con : q_constraints) {
        size_t n_nonzero_linear = con.lhs.l_coefs.size(),
                n_nonzero_quad = con.lhs.q_coefs.size();
        CPXDIM* linind = reinterpret_cast<CPXDIM*>(malloc(n_nonzero_linear * sizeof(CPXDIM)));
        double* linval = reinterpret_cast<double*>(malloc(n_nonzero_linear * sizeof(CPXDIM)));

        CPXDIM* quadrow = reinterpret_cast<CPXDIM*>(malloc(n_nonzero_quad * sizeof(CPXDIM)));
        CPXDIM* quadcol = reinterpret_cast<CPXDIM*>(malloc(n_nonzero_quad * sizeof(CPXDIM)));
        double* quadval = reinterpret_cast<double*>(malloc(n_nonzero_quad * sizeof(double)));

        size_t k = 0;
        for (auto& kv : con.lhs.l_coefs) {
            linind[k] = kv.first.column;
            linval[k] = kv.second;
            k++;
        }
        k = 0;
        for (auto& kv : con.lhs.q_coefs) {
            quadrow[k] = kv.first.v1.column;
            quadrow[k] = kv.first.v2.column;
            quadval[k] = kv.second;
            k++;
        }

        status = CPXXaddqconstr(env, prog, n_nonzero_linear, n_nonzero_quad, con.rhs, get_sense(con),
                                linind, linval, quadrow, quadcol, quadval, NULL);
        handle_status(status);
    }
}

template <class T> inline int
CPXLPManager<T>::solve(double* obj_p, int* solstat_p) {
    int status;
    if (prog_type == problem_type::lp) {
        status = CPXXmipopt(env, prog);
    } else if (prog_type == problem_type::qp || prog_type == problem_type::qcp ) {
//      status = CPXXqpopt(env, prog);
        status = CPXXmipopt(env, prog);
    } else if (prog_type == problem_type::mip) {
        status = CPXXmipopt(env, prog);
    } else {
        std::cerr << "CPXLPManager: unsupported problem type\n";
        return -1;
    }
    std::cout << "opt status: " << status << "\n";
    handle_status(status);

    int fincols = CPXXgetnumcols(env, prog);
    if (prog_soln != NULL) free(prog_soln);
    prog_soln = reinterpret_cast<double*>(malloc(fincols * sizeof(double)));

    status = CPXXsolution(env, prog, solstat_p, obj_p, prog_soln, NULL, NULL, NULL);
    return status;
}

template <class T> inline void
CPXLPManager<T>::solve_pool() {
    int status = CPXXpopulate(env, prog);
    handle_status(status);
}

template <class T> inline lp_var_t
CPXLPManager<T>::get_var(T label) {
    return label_to_lp_var[label];
}

template <class T> inline double
CPXLPManager<T>::get_value(T label) {
    return prog_soln[label_to_lp_var[label].column];
}

template <class T> inline lp_var_t
CPXLPManager<T>::add_slack_var(double lwr, double upp, lp_var_t::bounds btype, lp_var_t::domain vtype) {
    lp_var_t v;
    v.column = variables.size();
    v.lwr = lwr;
    v.upp = upp;
    v.bounds_type = btype;
    v.var_type = vtype;
    variables.push_back(v);
    return v;
}

template <class T> inline lp_var_t
CPXLPManager<T>::add_var(T label, double lwr, double upp, lp_var_t::bounds btype, lp_var_t::domain vtype) {
    lp_var_t v = add_slack_var(lwr, upp, btype, vtype);
    label_to_lp_var[label] = v;
    return v;
}

template <class T> inline void
CPXLPManager<T>::add_constraint(lp_constr_t con) {
    auto& constraints_array = con.is_quadratic() ? q_constraints : l_constraints;
    if (con.relation == lp_constr_t::direction::neq) {
        const double M = 10'000'000;
        // We will need to expand this constraint into multiple constraints.
        // To do so, we will need to introduce a slack variable y.
        lp_var_t y = add_slack_var(0, 1, lp_var_t::bounds::both, lp_var_t::domain::binary);
        // Want to add constraints:
        // lhs - My <= rhs - 1
        // lhs - My >= rhs + 1 - M
        constraints_array.emplace_back(con.lhs - M*y, con.rhs - 1, lp_constr_t::direction::le);
        constraints_array.emplace_back(con.lhs - M*y, con.rhs + 1 - M, lp_constr_t::direction::ge);
    } else {
        constraints_array.push_back(con);
    }
}

template <class T> inline size_t
CPXLPManager<T>::get_soln_pool_size() {
    return CPXXgetsolnpoolnumsolns(env, prog);
}

template <class T> inline double
CPXLPManager<T>::fetch_soln_from_pool(size_t k) {
    int fincols = CPXXgetnumcols(env, prog);
    if (prog_soln != NULL)  free(prog_soln);
    prog_soln = reinterpret_cast<double>(malloc(fincols * sizeof(double)));

    double obj;
    handle_status(CPXXgetsolnpoolobjval(env, prog, k, &obj));
    handle_status(CPXXgetsolnpoolx(env, prog, k, prog_soln, 0, fincols-1));
    return obj;
}

template <class T> inline char
CPXLPManager<T>::get_sense(lp_constr_t con) {
    if (con.relation == lp_constr_t::direction::ge) {
        return 'G';
    } else if (con.relation == lp_constr_t::direction::le) {
        return 'L';
    } else {
        return 'E';
    }
}

