/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

namespace vtils {

template <class SET_TYPE> inline SET_TYPE
immd_set_union(SET_TYPE s1, SET_TYPE s2) {
    for (const auto& x : s2) s1.insert(x);
    return s1;
}

template <class SET_TYPE> inline SET_TYPE
immd_set_diff(SET_TYPE s1, SET_TYPE s2) {
    for (const auto& x : s2) s1.erase(x);
    return s1;
}

template <class SET_TYPE> inline SET_TYPE
immd_set_intersect(SET_TYPE s1, SET_TYPE s2) {
    SET_TYPE s3;
    for (const auto& x : s1) {
        if (s2.count(x)) {
            s3.insert(x);
        }
    }
    return s3;
}

template <class SET_TYPE> inline SET_TYPE
immd_set_symdiff(SET_TYPE s1, SET_TYPE s2) {
    SET_TYPE s3;
    for (const auto& x : s1) {
        if (!s2.count(x)) {
            s3.insert(x);
        } else {
            s2.erase(x);
        }
    }
    upd_set_union(s3,s2);
    return s3;
}

template <class SET_TYPE> inline SET_TYPE&
upd_set_union(SET_TYPE& s1, SET_TYPE s2) { s1 = immd_set_union(s1,s2); return s1; }

template <class SET_TYPE> inline SET_TYPE&
upd_set_diff(SET_TYPE& s1, SET_TYPE s2) { s1 = immd_set_diff(s1,s2); return s1; }

template <class SET_TYPE> inline SET_TYPE&
upd_set_intersect(SET_TYPE& s1, SET_TYPE s2) { s1 = immd_set_intersect(s1,s2); return s1; }

template <class SET_TYPE> inline SET_TYPE&
upd_set_symdiff(SET_TYPE& s1, SET_TYPE s2) { s1 = immd_set_symdiff(s1,s2); return s1; }

}   // vtils
