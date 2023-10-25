/* author: Suhas Vittal
 *  date:   24 October 2023
 * */

#ifndef QONTRA_LINPROG_h
#define QONTRA_LINPROG_h

#include "defs.h"

#include <map>
#include <vector>

#include <cplex.h>

#include <stdlib.h>

namespace qontra {

template <typename T> class LPManager;

template <typename T> struct lp_var_t;
template <typename T> struct lp_expr_t;
template <typename T> struct lp_constr_t;

template <typename T>
struct lp_expr_t {
    lp_expr_t(void)
        :coefs(),
        constant(0)
    {}

    lp_expr_t(fp_t c)
        :coefs(),
        constant(c)
    {}

    lp_expr_t(lp_var_t<T>* v) 
        :coefs(),
        constant(0)
    {
        coefs[v] = 1.0;
    }

    lp_expr_t(const lp_expr_t<T>& other) 
        :coefs(other.coefs),
        constant(other.constant)
    {}

    lp_expr_t<T>& operator+=(const lp_expr_t<T> other) {
        constant += other.constant;
        for (auto& kv : coefs)  kv.second += other.coefs.at(kv.first);
        return *this;
    }

    lp_expr_t<T>& operator-=(const lp_expr_t<T> other) {
        constant -= other.constant;
        for (auto& kv : coefs)  kv.second -= other.coefs.at(kv.first);
        return *this;
    }

    lp_expr_t<T>& operator+=(fp_t x) {
        constant += x;
        return *this;
    }

    lp_expr_t<T>& operator-=(fp_t x) {
        constant -= x;
        return *this;
    }

    lp_expr_t<T>& operator*=(fp_t x) {
        for (auto& kv : coefs)  kv.second *= x;
        return *this;
    }

    lp_expr_t<T>& operator/=(fp_t x) {
        for (auto& kv : coefs)  kv.second /= x;
        return *this;
    }

    fp_t constant;
    std::map<lp_var_t<T>*, fp_t>   coefs;
};

template <typename T> lp_expr_t<T> operator+(lp_expr_t<T> e, lp_expr_t<T> f) { return (e += f); }
template <typename T> lp_expr_t<T> operator-(lp_expr_t<T> e, lp_expr_t<T> f) { return (e -= f); }

template <typename T> lp_expr_t<T> operator+(lp_expr_t<T> e, fp_t x) { return (e += x); }
template <typename T> lp_expr_t<T> operator-(lp_expr_t<T> e, fp_t x) { return (e -= x); }
template <typename T> lp_expr_t<T> operator*(lp_expr_t<T> e, fp_t x) { return (e *= x); }
template <typename T> lp_expr_t<T> operator/(lp_expr_t<T> e, fp_t x) { return (e /= x); }
template <typename T> lp_expr_t<T> operator+(fp_t x, lp_expr_t<T> e) { return e+x; }
template <typename T> lp_expr_t<T> operator-(fp_t x, lp_expr_t<T> e) { return e-x; }
template <typename T> lp_expr_t<T> operator*(fp_t x, lp_expr_t<T> e) { return e*x; }
template <typename T> lp_expr_t<T> operator/(fp_t x, lp_expr_t<T> e) { return e/x; }

enum class VarBounds { none, lower, upper, both, fixed };
enum class VarDomain { continuous, integer, binary };

template <typename T>
struct lp_var_t {
    LPManager<T>*   lp_ptr;

    int32_t     column;

    fp_t lwr;
    fp_t upp;

    VarBounds bounds_type;
    VarDomain var_type;
};

/*
template <typename T> lp_expr_t<T> operator+(lp_var_t<T>* v, lp_var_t<T>* w) { return lp_expr_t<T>(v) + lp_expr_t<T>(w); }
template <typename T> lp_expr_t<T> operator-(lp_var_t<T>* v, lp_var_t<T>* w) { return lp_expr_t<T>(v) - lp_expr_t<T>(w); }
template <typename T> lp_expr_t<T> operator+(lp_var_t<T>* v, fp_t x) { return lp_expr_t<T>(v) + x; }
template <typename T> lp_expr_t<T> operator-(lp_var_t<T>* v, fp_t x) { return lp_expr_t<T>(v) - x; }
template <typename T> lp_expr_t<T> operator*(lp_var_t<T>* v, fp_t x) { return lp_expr_t<T>(v) * x; }
template <typename T> lp_expr_t<T> operator/(lp_var_t<T>* v, fp_t x) { return lp_expr_t<T>(v) / x; }
template <typename T> lp_expr_t<T> operator+(fp_t x, lp_var_t<T>* v) { return v + x; }
template <typename T> lp_expr_t<T> operator-(fp_t x, lp_var_t<T>* v) { return v - x; }
template <typename T> lp_expr_t<T> operator*(fp_t x, lp_var_t<T>* v) { return v * x; }
template <typename T> lp_expr_t<T> operator/(fp_t x, lp_var_t<T>* v) { return v / x; }
*/

enum class ConstraintDirection { ge, le, eq, neq };

template <typename T>
struct lp_constr_t {
    lp_constr_t(lp_expr_t<T> _lhs, lp_expr_t<T> _rhs, ConstraintDirection r) 
        :lhs(_lhs),
        rhs(0),
        relation(r)
    {
        for (auto kv : _rhs.coefs) {
            lhs.coefs[kv.first] -= kv.second;
        }
        rhs = _rhs.constant - _lhs.constant;
        lhs.constant = 0;
    }

    lp_constr_t(lp_expr_t<T> _lhs, fp_t _rhs) 
        :lhs(_lhs),
        rhs(_rhs)
    {
        rhs -= _lhs.constant;
        lhs.constant = 0;
    }

    lp_expr_t<T> lhs;
    fp_t rhs;
    ConstraintDirection relation;
};

inline CPXENVptr cpxinit(void) {
    int status;
    CPXENVptr env = CPXopenCPLEX(&status);
    return env;
}

inline void cpxexit(CPXENVptr* env_p) {
    CPXcloseCPLEX(env_p);
}

template <typename T>
class LPManager {
public:
    LPManager(CPXENVptr _env=NULL)
        :variables(),
        constraints(),
        label_to_lp_var(),
        columns(0),
        rows(0),
        env(_env),
        prog(NULL),
        prog_soln(NULL),
        env_is_initialized_by_object(false),
        status(0)
    {
        env_is_initialized_by_object = (env==NULL);
        if (env == NULL) {
            env = cpxinit();
        }
    }

    ~LPManager() {
        for (auto v : variables) delete v;
        cpxfreeprog();
        if (env_is_initialized_by_object) {
            cpxexit(&env);
        }
    }

    void build(lp_expr_t<T> objective, bool is_maximization=true) {
        cpxfreeprog();
        cpxmakeprog();

        status = CPXchgobjsen(env, prog, is_maximization ? CPX_MAX : CPX_MIN);
        
        // Declare variables in columns.
        fp_t* obj = (fp_t*) calloc(columns, sizeof(fp_t));
        fp_t* lb = (fp_t*) calloc(columns, sizeof(fp_t));
        fp_t* ub = (fp_t*) calloc(columns, sizeof(fp_t));
        char* vtypes = (char*) calloc(columns, sizeof(char));

        for (int32_t i = 0; i < columns; i++) {
            lp_var_t<T>* v = variables[i];
            obj[i] = objective.coefs[v];

            bool lwr_defined = v->bounds_type == VarBounds::lower
                                || v->bounds_type == VarBounds::both
                                || v->bounds_type == VarBounds::fixed;
            bool upp_defined = v->bounds_type == VarBounds::upper
                                || v->bounds_type == VarBounds::both
                                || v->bounds_type == VarBounds::fixed;
            
            lb[i] = lwr_defined ? v->lwr : -CPX_INFBOUND;
            ub[i] = upp_defined ? v->upp : CPX_INFBOUND;
            if (v->var_type == VarDomain::continuous) {
                vtypes[i] = 'C';
            } else if (v->var_type == VarDomain::integer) {
                vtypes[i] = 'I';
            } else {
                vtypes[i] = 'B';
            }
        }
        status = CPXnewcols(env, prog, columns, obj, lb, ub, vtypes, NULL);
        // Declare constraints in rows.
        //
        // First, we need to count the total number of nonzero coefficients in 
        // the constraints.
        int num_nonzeros = 0;
        for (auto& con : constraints) {
            num_nonzeros += con.lhs.coefs.size();
        }
        // Now, we need to populate the data structures to declare the
        // constraints.
        fp_t* rhs = (fp_t*) malloc(rows * sizeof(fp_t));
        char* sense = (char*) malloc(rows * sizeof(char));
        int* rmatbeg = (int*) malloc(rows * sizeof(int));
        int* rmatind = (int*) malloc(num_nonzeros * sizeof(int));
        fp_t* rmatval = (fp_t*) malloc(num_nonzeros * sizeof(fp_t));

        int offset = 0;
        for (int32_t i = 0; i < rows; i++) {
            lp_constr_t<T> con = constraints[i];

            rhs[i] = con.rhs;
            if (con.relation == ConstraintDirection::ge) {
                sense[i] = 'G';
            } else if (con.relation == ConstraintDirection::le) {
                sense[i] = 'L';
            } else {
                sense[i] = 'E';
            }
            // We will not have anything else.

            rmatbeg[i] = offset;
            for (auto kv : con.lhs.coefs) {
                lp_var_t<T>* v = kv.first;
                rmatind[offset] = v->column;
                rmatval[offset] = kv.second;
                offset++;
            }
        }
        status = CPXaddrows(env, prog, 0, rows, num_nonzeros, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
        // Free all variables.
        free(obj);
        free(lb);
        free(ub);
        free(vtypes);
        free(rhs);
        free(sense);
        free(rmatbeg);
        free(rmatind);
        free(rmatval);
    }

    fp_t solve() {
        status = CPXlpopt(env, prog);

        int finrows = CPXgetnumrows(env, prog);
        int fincols = CPXgetnumcols(env, prog);

        int solstat;
        fp_t objective;

        if (prog_soln != NULL)  free(prog_soln);
        fp_t* prog_soln = (fp_t*) malloc(fincols * sizeof(fp_t));

        fp_t* slack = (fp_t*) malloc(finrows * sizeof(fp_t));
        fp_t* dj = (fp_t*) malloc(fincols * sizeof(fp_t));
        fp_t* pi = (fp_t*) malloc(finrows * sizeof(fp_t));

        status = CPXsolution(env, prog, &solstat, &objective, prog_soln, pi, slack, dj);
        free(pi);
        free(slack);
        free(dj);

        return objective;
    }

    lp_var_t<T>* get_var(T label) {
        return label_to_lp_var[label];
    }

    fp_t get_value(T label) {
        int c = label_to_lp_var[label]->column;
        return prog_soln[c];
    }

    lp_var_t<T>* add_slack_var(fp_t lower, fp_t upper, VarBounds btype, VarDomain vtype) {
        lp_var_t<T>* v = new lp_var_t<T>;
        v->lp_ptr = this;
        v->column = columns++;
        v->lwr = lower;
        v->upp = upper;
        v->bounds_type = btype;
        v->var_type = vtype;
        return v;
    }

    lp_var_t<T>* add_var(T label, fp_t lower, fp_t upper, VarBounds btype, VarDomain vtype) {
        auto v = add_slack_var(lower, upper, btype, vtype);
        label_to_lp_var[label] = v;
        return v;
    }

    int32_t add_constraint(lp_constr_t<T> con) {
        if (con.relation == ConstraintDirection::neq) {
            const fp_t M = 10'000'000;
            // We will need to expand this constraint into multiple constraints.
            // To do so, we will need to introduce a slack variable y.
            lp_var_t<T>* y = add_slack_var(0, 1, VarBounds::both, VarDomain::binary);
            // Want to add constraints:
            // lhs - My <= rhs-1
            // lhs - My >= rhs+1-M
            lp_constr_t<T> con1(con.lhs - M*lp_expr_t<T>(y), con.rhs-1, ConstraintDirection::le);
            lp_constr_t<T> con2(con.lhs - M*lp_expr_t<T>(y), con.rhs+1-M, ConstraintDirection::ge);
            constraints.push_back(con1);
            constraints.push_back(con2);
            return (rows += 2);
        } else {
            constraints.push_back(con);
            return (++rows);
        }
    }

    std::vector<lp_var_t<T>*>       variables;
    std::vector<lp_constr_t<T>>     constraints;
private:
    void cpxfreeprog(void) {
        if (prog != NULL) {
            status = CPXfreeprob(env, &prog);
            prog = NULL;
        }
    }

    void cpxmakeprog(void) {
        prog = CPXcreateprob(env, &status, "lp");
    }

    std::map<T, lp_var_t<T>*> label_to_lp_var;

    int32_t columns;
    int32_t rows;

    CPXENVptr env;
    CPXLPptr prog;
    bool env_is_initialized_by_object;
    fp_t* prog_soln;

    int status;
};

}   // qontra

#endif  // QONTRA_LINPROG_h
