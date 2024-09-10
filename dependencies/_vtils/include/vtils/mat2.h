/*
 *  author: Suhas Vittal
 *  date:   29 June 2024
 * */

#ifndef VTILS_MAT2_h
#define VTILS_MAT2_h

#include <iostream>
#include <vector>

#include <stdint.h>

// Matrix operations in Z_2

namespace vtils {

class Mat2 {
public:
    Mat2(size_t, size_t);
    Mat2(const Mat2&);
    ~Mat2(void);

    static Mat2 from_rows(std::vector<Mat2>);

    bool get(size_t, size_t) const;
    void set(size_t, size_t, bool);

    void set_row(size_t, Mat2);
    void copy_rows_from(size_t ri, Mat2 from);
    void copy_rows_from(size_t ri, std::vector<Mat2>);

    Mat2 get_row(size_t) const;
    Mat2 get_col(size_t) const;

    Mat2 rref(std::vector<size_t>& pivots) const;
    size_t popcount(size_t row) const;
    Mat2 append_row(Mat2) const;

    Mat2 tr(void) const;
    const Mat2& trr(void) const;

    bool operator()(size_t, size_t) const;

    const size_t n_rows;
    const size_t n_cols;
private:
    Mat2(size_t, size_t, Mat2* trmat_p);

    void force_update_trmat(void); // Should be used when data is modified directly.

    const size_t n_words;

    std::vector<std::vector<uint64_t>> data;
    Mat2* trmat_p; // Transposed matrix.
    bool trmat_is_heapalloc;

    friend Mat2 operator+(Mat2, Mat2);
    friend Mat2 operator*(Mat2, Mat2);
    friend Mat2 operator&(Mat2, Mat2);
    friend Mat2 _rref(const Mat2&, std::vector<size_t>&);
};

std::ostream& operator<<(std::ostream&, Mat2);

Mat2 operator+(Mat2, Mat2);
Mat2 operator*(Mat2, Mat2);
Mat2 operator&(Mat2, Mat2);

std::vector<Mat2>   get_basis_vectors(const Mat2&);
std::vector<Mat2>   get_null_basis_vectors(const Mat2&);
size_t              rk(const Mat2&);


// As the actual rref operation in Mat2 requires force_update_trmat (which is expensive),
// we instead have a helper (_rref) that does the RREF. If we call rref itself, then we
// perform the force_update_trmat. However, operations like get_basis_vectors do not need
// the output's transpose.
Mat2 _rref(const Mat2&, std::vector<size_t>&);

}   // vtils

#endif
