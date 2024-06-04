/*
 *  author: Suhas Vittal
 *  date:   29 May 2024
 * */

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

}   // placc
