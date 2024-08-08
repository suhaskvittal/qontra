/*
 *  author: Suhas Vittal
 *  date:   7 August 2024
 * */

namespace cgen {

inline void
link(sptr<check_t> c1, int c1s, sptr<check_t> c2, int c2s) {
    if (c1->get(c1s) != nullptr || c2->get(c2s) != nullptr) {
        std::cerr << "Writing to already allocated side.\n";
        exit(1);
    }
    c1->get(c1s) = c2;
    c2->get(c2s) = c1;
    c1->side_map[c2][c2s] = c1s;
    c2->side_map[c1][c1s] = c2s;
}

inline sptr<check_t>&
Tiling::operator()(int i, int j) {
    return in_plane[i][j];
}

inline sptr<check_t>
Tiling::at(int i, int j) const {
    return in_plane.at(i).at(j);
}

inline sptr<check_t>
Tiling::add_check_at(int color, int sides, int i, int j) {
    sptr<check_t> c = std::make_shared<check_t>(color, sides, i, j);
    if (i < 0 || j < 0) {
        out_of_plane.push_back(c);
    } else {
        in_plane[i][j] = c;
    }
    all.push_back(c);
    return c;
}

inline std::vector<sptr<check_t>>
Tiling::get_all_checks() const {
    return all;
}

}   // cgen
