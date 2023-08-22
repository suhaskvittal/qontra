/*
 *  author: Suhas Vittal
 * */

#include "defs.h"

#include <Bond.hpp>
#include <UniTensor.hpp>

#include <iostream>

#include <math.h>

using namespace qontra;

int main() {
    // GHZ test.

    cytnx::Tensor q_base({1, 2}, cytnx::Type.ComplexDouble);
    q_base.at({0, 0}) = 1;

    cytnx::Tensor h_base({2, 2}, cytnx::Type.ComplexDouble);
    h_base.at({0, 0}) = 1;
    h_base.at({0, 1}) = 1;
    h_base.at({1, 0}) = 1;
    h_base.at({1, 1}) = -1;
    h_base = h_base * sqrt(0.5);

    cytnx::Tensor cx_base({2, 2, 2, 2}, cytnx::Type.ComplexDouble);
    cx_base.at({0, 0, 0, 0}) = 1;
    cx_base.at({0, 1, 0, 1}) = 1;
    cx_base.at({1, 0, 1, 1}) = 1;
    cx_base.at({1, 1, 1, 0}) = 1;

    cytnx::UniTensor q1(q_base);
    q1.set_labels({std::string("-1"), std::string("1")});
    q1.set_name("q1");

    cytnx::UniTensor q2(q_base);
    q2.set_labels({std::string("-2"), std::string("2")});
    q2.set_name("q2");

    cytnx::UniTensor q3(q_base);
    q3.set_labels({std::string("-3"), std::string("3")});
    q3.set_name("q3");

    cytnx::UniTensor h1(h_base);
    h1.set_labels({std::string("1"), std::string("4")});
    h1.set_name("h");

    cytnx::UniTensor cx1(cx_base);
    cx1.set_labels({std::string("4"),
                        std::string("2"), 
                        std::string("5"), 
                        std::string("6")});
    cx1.set_name("cx");
    
    cytnx::UniTensor cx2(cx_base);
    cx2.set_labels({std::string("6"),
                        std::string("3"), 
                        std::string("7"), 
                        std::string("8")});
    cx2.set_name("cx");

    cytnx::UniTensor h2(h_base);
    h2.set_labels({std::string("5"), std::string("9")});

    auto s1 = cytnx::Contract(q1, h1);
    auto s2 = cytnx::Contract(s1, cx1);
    auto s3 = cytnx::Contract(q2, s2);
    auto s4 = cytnx::Contract(s3, cx2);
    auto s5 = cytnx::Contract(q3, s4);
    auto s6 = cytnx::Contract(s5, h2);
    s1.print_diagram();
    s2.print_diagram();
    s3.print_diagram();
    s4.print_diagram();
    s5.print_diagram();
    s6.print_diagram();
    s6.print_blocks();

    return 0;
}
