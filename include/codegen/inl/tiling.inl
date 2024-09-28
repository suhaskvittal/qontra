/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

namespace cgen {

inline sptr<check_t>&
Tiling::operator()(int i, int j) {
    return in_plane[i][j];
}

inline sptr<check_t>
Tiling::at(int i, int j) const {
    return in_plane.at(i).at(j);
}

inline std::vector<sptr<check_t>>
Tiling::get_checks(int c) const {
    return checks_by_color.at(c);
}

inline std::vector<sptr<check_t>>
Tiling::get_all_checks() const {
    return all;
}

inline const std::vector<sptr<check_t>>&
Tiling::get_checks_ref(int c) const {
    return checks_by_color.at(c);
}

inline const std::vector<sptr<check_t>>&
Tiling::get_all_checks_ref() const {
    return all;
}

}   // cgen
