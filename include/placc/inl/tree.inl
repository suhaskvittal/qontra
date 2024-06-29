/*
 *  author: Suhas Vittal
 *  date:   31 May 2024
 * */

namespace placc {
    
inline sptr<fpn_v_t>
ShorTree::get_head() const {
    return head;
}

inline std::vector<sptr<fpn_v_t>>
ShorTree::get_leaves() const {
    return leaves;
}

}   // placc
