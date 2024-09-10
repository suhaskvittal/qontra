/*
 *  author: Suhas Vittal
 *  date:   29 June 2024
 * */

#include "vtils/mat2.h"

#include <algorithm>

namespace vtils {

std::vector<std::vector<uint64_t>>
init_zeros(size_t n, size_t m) {
    std::vector<std::vector<uint64_t>> dat(
            n,
            std::vector<uint64_t>( ((m-1)>>6) + 1, 0 )
    );
    return dat;

}

Mat2::Mat2(size_t n, size_t m)
    :n_rows(n),
    n_cols(m),
    n_words( ((m-1)>>6)+1 ),
    data( init_zeros(n,m) ),
    trmat_p( new Mat2(m,n,this) ),
    trmat_is_heapalloc(true)
{}

Mat2::Mat2(const Mat2& m)
    :n_rows(m.n_rows),
    n_cols(m.n_cols),
    n_words(m.n_words),
    data(m.data),
    trmat_p( new Mat2(m.n_rows,m.n_cols,this) ),
    trmat_is_heapalloc(true)
{
    trmat_p->data = m.trmat_p->data;
}

Mat2::~Mat2() {
    if (trmat_is_heapalloc) delete trmat_p;
}

Mat2
Mat2::from_rows(std::vector<Mat2> rows) {
    Mat2 m( rows.size(), rows[0].n_cols );
    for (size_t i = 0; i < rows.size(); i++) {
        m.data[i] = rows[i].data[0];
    }
    return m;
}

bool
Mat2::get(size_t i, size_t j) const {
    return static_cast<bool>( 
                data.at(i).at(j>>6) & (1L << (j&63))
            );
}

void
Mat2::set(size_t i, size_t j, bool x) {
    if (x) {
        data[i][j>>6] |= 1L << (j&63);
        trmat_p->data[j][i>>6] |= 1L << (i&63);
    } else {
        data[i][j>>6] &= ~(1L << (j&63));
        trmat_p->data[j][i>>6] &= ~(1L << (i&63));
    }
}

void
Mat2::set_row(size_t i, Mat2 row) {
    data[i] = row.data[0];
    for (size_t j = 0; j < n_cols; j++) {
        trmat_p->set(j,i,get(i,j));
    }
}

void
Mat2::copy_rows_from(size_t ri, Mat2 from) {
    for (size_t i = 0; i < from.n_rows; i++) {
        data[i+ri] = from.data[i];
        for (size_t j = 0; j < n_cols; j++) {
            trmat_p->set(j,i+ri,get(i+ri,j));
        }
    }
}

void
Mat2::copy_rows_from(size_t ri, std::vector<Mat2> rows) {
    for (size_t i = 0; i < rows.size(); i++) {
        data[i+ri] = rows[i].data[0];
        for (size_t j = 0; j < n_cols; j++) {
            trmat_p->set(j,i+ri,get(i+ri,j));
        }
    }
}

Mat2
Mat2::get_row(size_t i) const {
    Mat2 r(1,n_cols);
    r.data[0] = data[i];
    return r;
}

Mat2
Mat2::get_col(size_t j) const {
    return trmat_p->get_row(j);
}

Mat2
Mat2::rref(std::vector<size_t>& pivots) const {
    Mat2 rr = _rref(*this, pivots);
    rr.force_update_trmat();
    return rr;
}

size_t
Mat2::popcount(size_t r) const {
    size_t c = 0;
    for (uint64_t w : data.at(r)) c += __builtin_popcountll(w);
    return c;
}

Mat2
Mat2::tr(void) const {
    return *trmat_p;
}

const Mat2&
Mat2::trr(void) const {
    return *trmat_p;
}

bool
Mat2::operator()(size_t i, size_t j) const {
    return get(i,j);
}

// PRIVATE FUNCTIONS

Mat2::Mat2(size_t n, size_t m, Mat2* t_p)
    :n_rows(n),
    n_cols(m),
    n_words( ((m-1)>>6) + 1 ),
    data( init_zeros(n,m) ),
    trmat_p(t_p),
    trmat_is_heapalloc(false)
{}

void
Mat2::force_update_trmat() {
    for (size_t j = 0; j < n_cols; j++) {
        for (size_t i = 0; i < n_rows; i++) {
            if (get(i,j)) {
                trmat_p->data[j][i>>6] |= 1L << (i&63);
            } else {
                trmat_p->data[j][i>>6] &= ~(1L << (i&63));
            }
        }
    }
}

// OPERATIONS

std::ostream&
operator<<(std::ostream& out, Mat2 m) {
    for (size_t i = 0; i < m.n_rows; i++) {
        if (i == 0) {
            out << "$[";
        } else {
            out << " [";
        }
        for (size_t j = 0; j < m.n_cols; j++) {
            if (m(i,j)) out << "!";
            else        out << ".";
        }
        if (i == m.n_rows-1) {
            out << "];" << std::endl;
        } else {
            out << "]" << std::endl;
        }
    }
    return out;
}

Mat2
operator+(Mat2 x, Mat2 y) {
    if (x.n_rows != y.n_rows && x.n_cols != y.n_cols) {
        std::cout << "[ Mat2 ] incompatible dimensions while adding" << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < x.n_rows; i++) {
        for (size_t j = 0; j < x.n_words; j++) {
            x.data[i][j] ^= y.data[i][j];
        }
    }
    x.force_update_trmat();
    return x;
}

Mat2
operator*(Mat2 x, Mat2 y) {
    if (x.n_cols != y.n_rows) {
        std::cout << "[ Mat2 ] incompatible dimensions while multiplying" << std::endl;
        exit(1);
    }

    Mat2 z(x.n_rows, y.n_cols);

    for (size_t i = 0; i < x.n_rows; i++) {
        for (size_t j = 0; j < y.n_cols; j++) {
            for (size_t k = 0; k < x.n_words; k++) {
                z.set(i,j,
                    static_cast<bool>(
                        __builtin_popcountll(x.data[i][k] & y.trr().data[j][k]) & 1
                        )
                    );
            }
        }
    }
    return z;
}

Mat2
operator&(Mat2 x, Mat2 y) {
    if (x.n_rows != y.n_rows && x.n_cols != y.n_cols) {
        std::cout << "[ Mat2 ] incompatible dimensions for elementwise multiply" << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < x.n_rows; i++) {
        for (size_t j = 0; j < x.n_words; j++) {
            x.data[i][j] &= y.data[i][j];
        }
    }
    x.force_update_trmat();
    return x;
}

std::vector<Mat2>
get_basis_vectors(const Mat2& m) {
    std::vector<size_t> p;
    _rref(m,p);
    std::vector<Mat2> vec;
    for (size_t i : p) {
        vec.push_back( m.get_col(i) );
    }
    return vec;
}

std::vector<Mat2>
get_null_basis_vectors(const Mat2& m) {
    std::vector<size_t> p;
    auto r = _rref(m,p);

    std::vector<Mat2> vec(m.n_cols - p.size(), Mat2(1, m.n_cols));
    size_t k = 0;
    for (size_t j = 0; j < m.n_cols; j++) {
        if (std::find(p.begin(), p.end(), j) != p.end()) continue;
        vec[k++].set(0, j, true);
    }
    for (size_t i = 0; i < p.size(); i++) {
        size_t jp = p.at(i);
        k = 0;
        for (size_t j = 0; j < m.n_cols; j++) {
            if (std::find(p.begin(), p.end(), j) != p.end()) continue;
            vec[k++].set(0, jp, r(i,j));
        }
    }
    return vec;
}

size_t
rk(const Mat2& m) {
    std::vector<size_t> p;
    _rref(m,p);
    return p.size();
}

Mat2
_rref(const Mat2& m, std::vector<size_t>& pivots) {
    Mat2 rr(m);

    size_t r = 0;
    for (size_t p = 0; p < m.n_cols && r < m.n_rows; p++) {
        if (!rr(r,p)) {
            size_t i = r+1;
            while ( i < m.n_rows && !rr(i,p) ) i++;
            if (i == m.n_rows) continue;
            std::swap( rr.data[r], rr.data[i] );
        }
        pivots.push_back(p);

        for (size_t i = r+1; i < m.n_rows; i++) {
            if (rr(i,p)) {
                for (size_t j = 0; j < m.n_words; j++) {
                    rr.data[i][j] ^= rr.data[r][j];
                }
            }
        }
        r++;
    }
    // rr is now in REF.
    // Now we must make rr become RREF.
    for (int i = pivots.size()-1; i >= 0; i--) {
        size_t j = pivots.at(i);
        for (size_t ii = 0; ii < i; ii++) {
            if (rr(ii,j)) {
                for (size_t jj = 0; jj < m.n_words; jj++) {
                    rr.data[ii][jj] ^= rr.data[i][jj];
                }
            }
        }
    }
    return rr;
}

}   // vtils
