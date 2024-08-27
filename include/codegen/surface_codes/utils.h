/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#ifndef CGEN_SC_UTILS_h
#define CGEN_SC_UTILS_h

#include <qontra/defs.h>

#include <array>

namespace cgen {

template <int R, int S>                             struct Qubit;
template <int R, int S, class TYPE, int MAX_SIZE>   class Check;
template <int R, int S>                             class Star;
template <int R, int S>                             class Face;

// Using raw pointers instead of shared for less overheads. Every
// bit of performance matters for this application.

template <int R, int S>
struct Qubit {
    Qubit()
        :x_ptr(0),
        z_ptr(0),
        x_checks(),
        z_checks()
    {}

    std::array<Star<R,S>*, 2> x_checks;
    std::array<Face<R,S>*, 2> z_checks;
    uint8_t x_ptr;
    uint8_t z_ptr;
};

template <int R, int S, class TYPE, int MAX_SIZE>
class Check {
public:
    bool add_qubit(sptr<Qubit<R,S>>);

    inline sptr<Qubit<R,S>>& operator[](size_t i) { return qubits[ isan(i) ]; }
    inline const sptr<Qubit<R,S>>& at(size_t i) const { return qubits.at( isan(i) ); }

    inline uint8_t get_size(void) { return size; }
    inline bool is_full(void) { return size == MAX_SIZE; }
protected:
    constexpr Check() 
        :qubits(), size(0)
    {
        qubits.fill(nullptr);
    }

    inline uint8_t isan(uint8_t i) const { return i % MAX_SIZE; }

    std::array<sptr<Qubit<R,S>>, MAX_SIZE> qubits;
    uint8_t size;
private:
    inline constexpr TYPE* _cast(void) { return static_cast<TYPE*>(this); }
};

template <int R, int S>
class Star : public Check<R, S, Star<R,S>, R> {
public:
    Star()
        :Check<R,S,Star<R,S>,R>()
    {
        for (uint8_t i = 0; i < R; i++) this->add_qubit( std::make_shared<Qubit<R,S>>() );
    }

    void tie(Star<R,S>* other, uint8_t i, uint8_t j);
    
    bool update_qubit(sptr<Qubit<R,S>>);
private:
};

template <int R, int S>
class Face : public Check<R, S, Face<R,S>, S> {
public:
private:
    bool update_qubit(sptr<Qubit<R,S>>) const;
};

}   // cgen

#include "utils.inl"

#endif  // CGEN_SC_UTILS_h
