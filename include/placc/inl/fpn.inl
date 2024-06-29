/*
 *  author: Suhas Vittal
 *  date:   29 May 2024
 * */

namespace qontra {
namespace graph {

template <> inline std::string
print_v<placc::fpn_v_t>(sptr<placc::fpn_v_t> v) {
    using namespace placc;

    if (v == nullptr) return "<nil>";

    std::string s;
    if (v->qubit_type == fpn_v_t::type::data) {
        s += "d";
    } else if (v->qubit_type == fpn_v_t::type::parity) {
        s += "p";
    } else if (v->qubit_type == fpn_v_t::type::flag) {
        s += "f";
    }
    s += std::to_string(v->id);
    return s;
}

}   // graph
}   // qontra

namespace placc {

inline void
push_back_gate(
        qes::Program<>& prog,
        std::string name,
        const std::vector<sptr<fpn_v_t>>& operands) 
{
    std::vector<int64_t> _operands;
    for (sptr<fpn_v_t> x : operands) _operands.push_back( static_cast<int64_t>(x->id) );
    prog.emplace_back(name, _operands);
}

inline void
push_back_gate(
        qes::Program<>& prog,
        std::string name,
        std::initializer_list<std::vector<sptr<fpn_v_t>>> operand_list)
{
    std::vector<int64_t> _operands;
    for (const std::vector<sptr<fpn_v_t>>& operands : operand_list) {
        for (sptr<fpn_v_t> x : operands) _operands.push_back( static_cast<int64_t>(x->id) );
    }
    prog.emplace_back(name, _operands);
}

inline void
push_back_measurement(
        qes::Program<>& prog,
        const std::vector<sptr<fpn_v_t>>& operands,
        size_t& ctr,
        mm_t& m)
{
    std::vector<int64_t> _operands;
    for (const sptr<fpn_v_t>& x : operands) {
        m[x] = ctr++;
        _operands.push_back( static_cast<int64_t>(x->id) );
    }
    prog.emplace_back("measure", _operands);
}

inline void
push_back_measurement(
        qes::Program<>& prog,
        std::initializer_list<std::vector<sptr<fpn_v_t>>> operand_list,
        size_t& ctr,
        mm_t& m)
{
    std::vector<int64_t> _operands;
    for (const std::vector<sptr<fpn_v_t>>& operands : operand_list) {
        for (const sptr<fpn_v_t>& x : operands) {
            m[x] = ctr++;
            _operands.push_back( static_cast<int64_t>(x->id) );
        }
    }
    prog.emplace_back("measure", _operands);
}

}   // placc
